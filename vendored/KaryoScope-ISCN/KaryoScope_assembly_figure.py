#!/usr/bin/env python3
# Vendored from barthel-lab/KaryoScope-ISCN at commit 4e42a98d6526ecb31ca6fb5d6c8910e4e84c605a, file KaryoScope_assembly_figure.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""
KaryoScope Assembly Figure - Combined letter-format visualization.

Combines all chromosomes (TACCAT genome plot) with featureset hierarchy
legends on a single portrait letter page.

Layout:
    Row 1: chr1-chr8
    Row 2: chr9-chr16
    Row 3: chr17-chr22, X+Y
    Bottom: chromosome / region / subtelomeric hierarchy trees

Usage:
    python KaryoScope_assembly_figure.py \\
        --bed-dir data/ \\
        --colors-dir /path/to/KS_human_CHM13/ \\
        --centromeres HG002.assembly.KS_human_CHM13.region.presmoothed.scaffolded.centromeres.bed \\
        --sample HG002.assembly \\
        --database KS_human_CHM13 \\
        --sex XY \\
        --output HG002_figure.svg
"""

import argparse
import subprocess
import sys
from pathlib import Path

import drawsvg as draw
from PIL import Image as PILImage

PROJECT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(PROJECT_DIR))
import KaryoScope_assembly_genome_plot as gp
import plot_featureset_tree_bars as tp

# plot_featureset_tree_bars.py's built-in defaults point to sibling directories
# in KaryoScope-ISCN (featuresets-hierarchy/, databases/KS_human_CHM13_v2/, ...)
# which is where this script also lives -- so no path overrides are needed.
# Users can override via --hierarchy-dir / --v2-db-dir if they have custom data.


def draw_scale_legend(d, x, y, pixels_per_bp_genome, cen_ppb_fixed,
                      telo_ppb_fixed, text_color):
    """Draw horizontal genome, centromere, and telomere scale bars + zoom ratios.

    Layout (top-left):
        Genome:     [====] 50 Mb
        Centromere: [====] 10 Mb  (4.9x zoom)
        Telomere:   [====] 25 kb  (2157x zoom)
    """
    label_font = 9
    bar_h = 2
    tick_h = 5
    label_x_pad = 8
    line_height = 18

    rows = [
        ("Genome:", 50_000_000, pixels_per_bp_genome, "Mb", False),
        ("Centromere:", 10_000_000, cen_ppb_fixed, "Mb", True),
        ("Telomere:", 25_000, telo_ppb_fixed, "kb", True),
    ]

    label_width = 70  # space reserved for the row labels
    bar_x = x + label_width

    for i, (row_label, length_bp, ppb, unit, show_zoom) in enumerate(rows):
        ry = y + i * line_height
        d.append(draw.Text(row_label, font_size=label_font, x=x, y=ry + 3,
                           fill=text_color, font_family="sans-serif",
                           text_anchor="start"))
        bar_w = length_bp * ppb
        d.append(draw.Rectangle(bar_x, ry, bar_w, bar_h, fill=text_color))
        d.append(draw.Line(bar_x, ry - tick_h / 2 + bar_h / 2,
                           bar_x, ry + tick_h / 2 + bar_h / 2,
                           stroke=text_color, stroke_width=0.8))
        d.append(draw.Line(bar_x + bar_w, ry - tick_h / 2 + bar_h / 2,
                           bar_x + bar_w, ry + tick_h / 2 + bar_h / 2,
                           stroke=text_color, stroke_width=0.8))
        if unit == "Mb":
            length_str = f"{length_bp // 1_000_000} Mb"
        else:
            length_str = f"{length_bp // 1_000} kb"
        if show_zoom and pixels_per_bp_genome > 0:
            zoom = ppb / pixels_per_bp_genome
            label = f"{length_str}  ({zoom:.0f}\u00d7 zoom)" if zoom >= 100 \
                    else f"{length_str}  ({zoom:.1f}\u00d7 zoom)"
        else:
            label = length_str
        d.append(draw.Text(label, font_size=label_font,
                           x=bar_x + bar_w + label_x_pad, y=ry + 3,
                           fill=text_color, font_family="sans-serif",
                           text_anchor="start"))

    return line_height * len(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combined KaryoScope assembly figure")
    # Required data paths
    parser.add_argument("--bed-dir", required=True,
                        help="Directory containing pre-binned BED feature files")
    parser.add_argument("--colors-dir", required=True,
                        help="Directory with color files")
    parser.add_argument("--centromeres", required=True,
                        help="Centromeres BED file path")
    parser.add_argument("--sample", required=True,
                        help="Sample / assembly name (e.g. HG002.assembly)")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output SVG file")
    # Visual options
    parser.add_argument("--background", default="white",
                        choices=["white", "black"])
    parser.add_argument("--sex", default="XY")
    parser.add_argument("--scale", type=float, default=1.4,
                        help="Pixels per Mb for genome track")
    parser.add_argument("--subtelo-span", type=int, default=25_000)
    parser.add_argument("--layout", default="grid", choices=["grid", "packed"],
                        help="Layout: grid (3 rows + trees) or packed (shelf bin-packing)")
    # Hierarchy tree data overrides (default: ISCN's built-in data)
    parser.add_argument("--hierarchy-dir", default=None,
                        help="Override hierarchy (v1) directory (default: ISCN built-in)")
    parser.add_argument("--v2-db-dir", default=None,
                        help="Override V2 hierarchy database directory")
    parser.add_argument("--v2-output-dir", default=None,
                        help="Override V2 hierarchy output directory")
    parser.add_argument("--hierarchy-version", default="v1", choices=["v1", "v2"],
                        help="Which hierarchy version to load for legends (default: v1)")
    # G-banding overlay
    parser.add_argument("--gbanding-dir", default=None,
                        help="Directory with per-chromosome G-banding PNG images")
    parser.add_argument("--gbanding-position", default="above",
                        choices=["above", "below"],
                        help="Place G-banding images above or below centromere panels")
    return parser.parse_args()


# G-banding image constants
GBAND_PAD = 4          # vertical padding between g-band image and centromere panel


def get_gbanding_path(gbanding_dir, chrom):
    """Return the path to a g-banding PNG for a chromosome, or None."""
    if gbanding_dir is None:
        return None
    d = Path(gbanding_dir)
    # Try chrXX.png (zero-padded) and chrX.png
    num = chrom.replace("chr", "")
    for name in [f"chr{num.zfill(2)}.png", f"chr{num}.png"]:
        p = d / name
        if p.exists():
            return str(p)
    return None


def get_gbanding_height(gbanding_path, display_width):
    """Return display height for a g-banding image at given width, preserving aspect ratio."""
    im = PILImage.open(gbanding_path)
    return display_width * im.height / im.width


def draw_genome_row(d, chromosomes, chrom_draw_haps, chrom_map, scaffold_lengths,
                    centromere_coords, all_chrom_features, colors,
                    pixels_per_bp_genome, genome_height, subtelo_h, cen_ppb_fixed,
                    max_cen_panel_h, subtelo_span, sex, background,
                    x_offset, y_offset, chrom_gap,
                    gbanding_dir=None, gbanding_position="above"):
    """Draw a row of chromosomes at given position. Returns row width."""
    text_color = "black" if background == "white" else "white"
    haplotypes = ["hap1", "hap2"]

    def _cen_on_right(chrom, hap):
        if sex == "male" and chrom == "chrX" and hap.startswith("hap2"):
            return True
        if sex == "male" and chrom == "chrY" and hap.startswith("hap1"):
            return False
        return hap.startswith("hap1")

    # Compute X positions for this row
    cx = x_offset
    chrom_positions = {}  # chrom -> x_start
    hap_positions = {}    # (chrom, hap) -> hx
    hap_label_info = {}   # (chrom, hap) -> (label_x, label_w)

    offset_extra = gp.TRACK_W + gp.SUBTELO_OFFSET_GAP

    drawable_in_row = [c for c in chromosomes if c in chrom_draw_haps]
    for idx, chrom in enumerate(drawable_in_row):
        draw_haps = chrom_draw_haps[chrom]
        chrom_positions[chrom] = cx

        hx = cx
        for i, hap in enumerate(draw_haps):
            is_h1 = _cen_on_right(chrom, hap)
            # Check if offset needed
            scaf = chrom_map.get(chrom, {}).get(hap)
            has_off = False
            if scaf and scaf in scaffold_lengths:
                scaf_len = scaffold_lengths[scaf]
                scaf_h = scaf_len * pixels_per_bp_genome
                cen_s, cen_e = centromere_coords.get(scaf, (None, None))
                this_cen_h = (cen_e - cen_s) * cen_ppb_fixed if cen_s is not None else 0
                eff_h = max(scaf_h, this_cen_h)
                has_off = (subtelo_h * 2 > eff_h)

            hap_positions[(chrom, hap)] = hx
            if has_off and is_h1:
                lx = hx - offset_extra
                lw = gp.hap_width() + offset_extra
            elif has_off and not is_h1:
                lx = hx
                lw = gp.hap_width() + offset_extra
                hx += offset_extra
            else:
                lx = hx
                lw = gp.hap_width()

            hap_label_info[(chrom, hap)] = (lx, lw)
            hx += gp.hap_width()
            if i < len(draw_haps) - 1:
                hx += gp.HAP_GAP

        # No chromosome gap between chrX and chrY (treat them as a unit)
        next_chrom = drawable_in_row[idx + 1] if idx + 1 < len(drawable_in_row) else None
        if {chrom, next_chrom} == {"chrX", "chrY"}:
            cx = hx + gp.HAP_GAP
        else:
            cx = hx + chrom_gap

    row_width = cx - chrom_gap - x_offset

    # Header labels
    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        draw_haps = chrom_draw_haps[chrom]
        first_lx, first_lw = hap_label_info[(chrom, draw_haps[0])]
        last_lx, last_lw = hap_label_info[(chrom, draw_haps[-1])]
        cl_x = first_lx
        cl_w = (last_lx + last_lw) - first_lx

        d.append(draw.Rectangle(cl_x, y_offset + gp.CHROM_LINE_Y, cl_w, 1,
                                fill=text_color))
        d.append(draw.Text(chrom, font_size=14,
                           x=cl_x + cl_w / 2, y=y_offset + gp.CHROM_LABEL_Y,
                           fill=text_color, font_family="sans-serif",
                           text_anchor="middle"))

        for hap in draw_haps:
            lx, lw = hap_label_info[(chrom, hap)]
            d.append(draw.Rectangle(lx, y_offset + gp.HAP_LINE_Y, lw, 1,
                                    fill=text_color))
            hap_text = f"h{hap[3:]}" if hap.startswith("hap") else hap
            d.append(draw.Text(hap_text, font_size=14,
                               x=lx + lw / 2, y=y_offset + gp.HAP_LABEL_Y,
                               fill=text_color, font_family="sans-serif",
                               text_anchor="middle"))

    y0 = y_offset + gp.INITIAL_Y

    # Helper to compute positions (same logic as genome plot)
    def compute_positions(scaf, scaf_len, is_h1, hx):
        scaf_h = scaf_len * pixels_per_bp_genome
        if is_h1:
            base_stx = hx
            genome_x = hx + gp.TRACK_W + gp.FAN_GAP
            cen_x = genome_x + gp.GENOME_TRACK_W + gp.FAN_GAP
        else:
            cen_x = hx
            genome_x = hx + gp.TRACK_W + gp.FAN_GAP
            base_stx = genome_x + gp.GENOME_TRACK_W + gp.FAN_GAP

        cen_s, cen_e = centromere_coords.get(scaf, (None, None))
        this_cen_h = (cen_e - cen_s) * cen_ppb_fixed if cen_s is not None else 0

        if this_cen_h > scaf_h:
            gy_off = (this_cen_h - scaf_h) / 2
        else:
            gy_off = 0
        gy0 = y0 + gy_off

        cpy = y0
        if cen_s is not None:
            cen_mid = gy0 + ((cen_s + cen_e) / 2) * pixels_per_bp_genome
            cpy = cen_mid - this_cen_h / 2
            cpy = max(y0, min(y0 + genome_height - this_cen_h, cpy))

        if this_cen_h > scaf_h and cen_s is not None:
            eff_top, eff_bot = cpy, cpy + this_cen_h
        else:
            eff_top, eff_bot = gy0, gy0 + scaf_h

        p_y = eff_top
        q_y = max(eff_top, eff_bot - subtelo_h)

        p_stx = base_stx
        if p_y + subtelo_h > q_y:
            if is_h1:
                p_stx = base_stx - gp.TRACK_W - gp.SUBTELO_OFFSET_GAP
            else:
                p_stx = base_stx + gp.TRACK_W + gp.SUBTELO_OFFSET_GAP

        return {
            "genome_y0": gy0, "p_panel_y": p_y, "q_panel_y": q_y,
            "p_stx": p_stx, "q_stx": base_stx,
            "genome_x": genome_x, "cen_x": cen_x,
            "cen_panel_y": cpy, "cen_panel_h": this_cen_h,
        }

    # Pass 1: fans
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
            hx = hap_positions[(chrom, hap)]
            is_h1 = _cen_on_right(chrom, hap)
            pos = compute_positions(scaf, scaf_len, is_h1, hx)
            gy0 = pos["genome_y0"]

            # p-arm fan
            gp_top = gy0
            gp_bot = gy0 + subtelo_span * pixels_per_bp_genome
            if gp_bot - gp_top < 2: gp_bot = gp_top + 2
            if is_h1:
                gp.draw_zoom_fan(d, gp_top, gp_bot, pos["p_panel_y"], pos["p_panel_y"] + subtelo_h,
                                 pos["genome_x"], pos["p_stx"] + gp.TRACK_W, direction="left")
            else:
                gp.draw_zoom_fan(d, gp_top, gp_bot, pos["p_panel_y"], pos["p_panel_y"] + subtelo_h,
                                 pos["genome_x"] + gp.GENOME_TRACK_W, pos["p_stx"], direction="right")

            # q-arm fan
            q_start = max(0, scaf_len - subtelo_span)
            gq_top = gy0 + q_start * pixels_per_bp_genome
            gq_bot = gy0 + scaf_h
            if gq_bot - gq_top < 2: gq_top = gq_bot - 2
            if is_h1:
                gp.draw_zoom_fan(d, gq_top, gq_bot, pos["q_panel_y"], pos["q_panel_y"] + subtelo_h,
                                 pos["genome_x"], pos["q_stx"] + gp.TRACK_W, direction="left")
            else:
                gp.draw_zoom_fan(d, gq_top, gq_bot, pos["q_panel_y"], pos["q_panel_y"] + subtelo_h,
                                 pos["genome_x"] + gp.GENOME_TRACK_W, pos["q_stx"], direction="right")

            # centromere fan
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = pos["cen_panel_h"]
            if cen_s is not None and this_cen_h > 0:
                gc_top = gy0 + cen_s * pixels_per_bp_genome
                gc_bot = min(gy0 + cen_e * pixels_per_bp_genome, gy0 + scaf_h)
                cpy = pos["cen_panel_y"]
                if is_h1:
                    gp.draw_zoom_fan(d, gc_top, gc_bot, cpy, cpy + this_cen_h,
                                     pos["genome_x"] + gp.GENOME_TRACK_W, pos["cen_x"], direction="right")
                else:
                    gp.draw_zoom_fan(d, gc_top, gc_bot, cpy, cpy + this_cen_h,
                                     pos["genome_x"], pos["cen_x"] + gp.TRACK_W, direction="left")

    # Pass 2: tracks + labels
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
            hx = hap_positions[(chrom, hap)]
            is_h1 = _cen_on_right(chrom, hap)
            pos = compute_positions(scaf, scaf_len, is_h1, hx)
            gy0 = pos["genome_y0"]

            gp.draw_feature_track(d, feats["genome"], colors["chromosome"],
                                  pos["genome_x"], gy0, gp.GENOME_TRACK_W, 0, pixels_per_bp_genome)
            gp.draw_track_outline(d, pos["genome_x"], gy0, gp.GENOME_TRACK_W, scaf_h, text_color)

            p_ppb = subtelo_h / subtelo_span
            gp.draw_feature_track(d, feats["p_arm"], colors["subtelomeric"],
                                  pos["p_stx"], pos["p_panel_y"], gp.TRACK_W, 0, p_ppb)
            gp.draw_track_outline(d, pos["p_stx"], pos["p_panel_y"], gp.TRACK_W, subtelo_h, text_color)

            q_start = max(0, scaf_len - subtelo_span)
            q_ppb = subtelo_h / max(1, scaf_len - q_start)
            gp.draw_feature_track(d, feats["q_arm"], colors["subtelomeric"],
                                  pos["q_stx"], pos["q_panel_y"], gp.TRACK_W, q_start, q_ppb)
            gp.draw_track_outline(d, pos["q_stx"], pos["q_panel_y"], gp.TRACK_W, subtelo_h, text_color)

            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = pos["cen_panel_h"]
            if cen_s is not None and feats["centromere"] and this_cen_h > 0:
                cen_span = cen_e - cen_s
                cpy = pos["cen_panel_y"]
                cen_ppb = this_cen_h / cen_span
                gp.draw_feature_track(d, feats["centromere"], colors["region"],
                                      pos["cen_x"], cpy, gp.TRACK_W, cen_s, cen_ppb)
                gp.draw_track_outline(d, pos["cen_x"], cpy, gp.TRACK_W, this_cen_h, text_color)

            # Size labels
            label_pad = 4
            genome_label = gp.format_size_label(scaf_len)
            glx = pos["genome_x"] + gp.GENOME_TRACK_W / 2
            gly = gy0 + scaf_h + label_pad
            d.append(draw.Text(genome_label, font_size=10, x=glx, y=gly,
                               fill=text_color, font_family="sans-serif",
                               text_anchor="end", dominant_baseline="middle",
                               transform=f"rotate(-90 {glx} {gly})"))
            if cen_s is not None and this_cen_h > 0:
                cen_span = cen_e - cen_s
                cen_label = gp.format_size_label(cen_span)
                clx = pos["cen_x"] + gp.TRACK_W / 2
                cly = pos["cen_panel_y"] + this_cen_h + label_pad
                d.append(draw.Text(cen_label, font_size=10, x=clx, y=cly,
                                   fill=text_color, font_family="sans-serif",
                                   text_anchor="end", dominant_baseline="middle",
                                   transform=f"rotate(-90 {clx} {cly})"))

    # Pass 3: g-banding images (above or below centromere panels)
    if gbanding_dir:
        for chrom in chromosomes:
            if chrom not in chrom_draw_haps:
                continue
            draw_haps = chrom_draw_haps[chrom]
            if len(draw_haps) != 2:
                continue
            gb_path = get_gbanding_path(gbanding_dir, chrom)
            if not gb_path:
                continue

            # Compute centromere positions for both haplotypes
            cen_positions = []
            for hap in draw_haps:
                scaf = chrom_map.get(chrom, {}).get(hap)
                if not scaf or scaf not in scaffold_lengths:
                    continue
                scaf_len = scaffold_lengths[scaf]
                hx = hap_positions[(chrom, hap)]
                is_h1 = _cen_on_right(chrom, hap)
                pos = compute_positions(scaf, scaf_len, is_h1, hx)
                cen_positions.append(pos)

            if len(cen_positions) != 2:
                continue

            # G-band image spans the two centromere bars
            cen_x_left = min(p["cen_x"] for p in cen_positions)
            cen_x_right = max(p["cen_x"] for p in cen_positions) + gp.TRACK_W
            gb_w = cen_x_right - cen_x_left
            gb_h = get_gbanding_height(gb_path, gb_w)

            # Get centromere panel vertical extent
            cen_top = min(p["cen_panel_y"] for p in cen_positions)
            cen_bot = max(p["cen_panel_y"] + p["cen_panel_h"] for p in cen_positions)

            if gbanding_position == "above":
                gb_y = cen_top - GBAND_PAD - gb_h
            else:  # below
                gb_y = cen_bot + GBAND_PAD

            d.append(draw.Image(cen_x_left, gb_y, gb_w, gb_h,
                                path=gb_path, embed=True))

    return row_width


def main():
    args = parse_args()
    # Data paths come directly from CLI args.
    paths = {
        "bed_dir": args.bed_dir,
        "colors_dir": args.colors_dir,
        "centromeres": args.centromeres,
        "sample": args.sample,
        "database": args.database,
    }

    # Optional hierarchy-tree path overrides (modifies the tp module at runtime).
    if args.hierarchy_dir:
        tp.HIERARCHY_DIR = Path(args.hierarchy_dir)
        tp.OUTPUT_DIR = tp.HIERARCHY_DIR
    if args.v2_db_dir:
        tp.V2_DB_DIR = Path(args.v2_db_dir)
    if args.v2_output_dir:
        tp.V2_OUTPUT_DIR = Path(args.v2_output_dir)

    pixels_per_bp_genome = args.scale / 1_000_000
    background = args.background
    text_color = "black" if background == "white" else "white"
    sex = gp.normalize_sex(args.sex)

    # Load genome data
    featuresets = ["chromosome", "subtelomeric", "region"]
    colors = gp.load_colors(paths["colors_dir"], paths["database"], featuresets,
                            background=background)
    centromere_coords = gp.load_centromeres(paths["centromeres"])

    genome_bed_base = f"{paths['sample']}.{paths['database']}.chromosome.presmoothed.scaffolded.binned1000000.features.bed"
    genome_bed = str(Path(paths["bed_dir"]) / genome_bed_base)
    if not Path(genome_bed).exists():
        genome_bed += ".gz"
    scaffold_lengths, chrom_map = gp.discover_scaffolds(genome_bed)

    haplotypes = sorted(set(h for haps in chrom_map.values() for h in haps)) or ["hap1", "hap2"]
    all_chroms = sorted(chrom_map.keys(), key=gp.chromosome_sort_key)

    chrom_draw_haps = {}
    for chrom in all_chroms:
        hap_data = chrom_map.get(chrom, {})
        expected = gp.get_expected_haps(chrom, sex, haplotypes)
        actual = [h for h in haplotypes if h in hap_data]
        # Always draw actual haplotypes (don't show empty synthetic-hap
        # placeholders inherited from other chromosomes). For sex chroms
        # only fall back to expected when no data exists at all.
        if actual:
            draw_haps = actual
        elif chrom in ("chrX", "chrY"):
            draw_haps = [h for h in haplotypes if h in expected]
        else:
            draw_haps = []
        if draw_haps or actual:
            chrom_draw_haps[chrom] = draw_haps

    max_genome_len = max(
        scaffold_lengths.get(chrom_map.get(c, {}).get(h, ""), 0)
        for c in chrom_draw_haps for h in chrom_draw_haps[c]
    )
    genome_height = max_genome_len * pixels_per_bp_genome
    subtelo_h = genome_height * gp.SUBTELO_FRAC

    max_cen_span = max((ce - cs for cs, ce in centromere_coords.values()), default=0)
    max_cen_panel_h = genome_height * gp.CEN_FRAC
    cen_ppb_fixed = max_cen_panel_h / max_cen_span if max_cen_span > 0 else 0

    # Load all chromosome features
    print("Loading features...")
    all_chrom_features = {}
    for chrom in all_chroms:
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
            hap_info = {"scaffold": scaf, "scaffold_length": scaf_len,
                        "cen_start": cen_s, "cen_end": cen_e}
            chrom_feats[hap] = gp.load_scaffold_features(
                hap_info, paths["bed_dir"], paths["sample"],
                paths["database"], "presmoothed", args.subtelo_span)
        all_chrom_features[chrom] = chrom_feats

    # Load hierarchy trees
    print("Loading hierarchy trees...")
    if args.hierarchy_version == "v2":
        tree_fs = tp.load_all_featuresets_v2()
    else:
        tree_fs = tp.load_all_featuresets()
    for root in tree_fs.values():
        tp._apply_color_overrides(root)

    tree_panels = {}
    for fs_name in ["chromosome", "region", "subtelomeric"]:
        if fs_name in tree_fs:
            tree_panels[fs_name] = tp.draw_panel(tree_fs[fs_name], fs_name)

    MARGIN = 30

    if args.layout == "grid":
        # Original grid layout: 3 rows of chromosomes + 1 row of hierarchies
        draw_grid_layout(
            d=None,  # will create drawing inside
            args=args, chrom_draw_haps=chrom_draw_haps, chrom_map=chrom_map,
            scaffold_lengths=scaffold_lengths, centromere_coords=centromere_coords,
            all_chrom_features=all_chrom_features, colors=colors,
            pixels_per_bp_genome=pixels_per_bp_genome, genome_height=genome_height,
            subtelo_h=subtelo_h, cen_ppb_fixed=cen_ppb_fixed,
            max_cen_panel_h=max_cen_panel_h, sex=sex, background=background,
            all_chroms=all_chroms, tree_panels=tree_panels,
        )
        return
    else:
        draw_packed_layout(
            args=args, chrom_draw_haps=chrom_draw_haps, chrom_map=chrom_map,
            scaffold_lengths=scaffold_lengths, centromere_coords=centromere_coords,
            all_chrom_features=all_chrom_features, colors=colors,
            pixels_per_bp_genome=pixels_per_bp_genome, genome_height=genome_height,
            subtelo_h=subtelo_h, cen_ppb_fixed=cen_ppb_fixed,
            max_cen_panel_h=max_cen_panel_h, sex=sex, background=background,
            all_chroms=all_chroms, tree_panels=tree_panels,
        )
        return


def draw_grid_layout(d, args, chrom_draw_haps, chrom_map, scaffold_lengths,
                     centromere_coords, all_chrom_features, colors,
                     pixels_per_bp_genome, genome_height, subtelo_h,
                     cen_ppb_fixed, max_cen_panel_h, sex, background,
                     all_chroms, tree_panels):
    pixels_per_bp_genome = args.scale / 1_000_000
    text_color = "black" if background == "white" else "white"

    # Layout: 3 rows of chromosomes + 1 row of hierarchies
    rows = [
        [c for c in all_chroms if c in chrom_draw_haps][:8],
        [c for c in all_chroms if c in chrom_draw_haps][8:16],
        [c for c in all_chroms if c in chrom_draw_haps][16:],
    ]

    CHROM_GAP_FIG = 30
    MARGIN = 30
    ROW_GAP = 15
    TREE_GAP = 20

    # Compute genome row heights (max scaffold height per row)
    row_heights = []
    for row_chroms in rows:
        max_h = 0
        for c in row_chroms:
            for h in chrom_draw_haps.get(c, []):
                scaf = chrom_map.get(c, {}).get(h)
                if scaf and scaf in scaffold_lengths:
                    scaf_len = scaffold_lengths[scaf]
                    scaf_h = scaf_len * pixels_per_bp_genome
                    cen_s, cen_e = centromere_coords.get(scaf, (None, None))
                    this_cen_h = (cen_e - cen_s) * cen_ppb_fixed if cen_s is not None else 0
                    max_h = max(max_h, scaf_h, this_cen_h)
        row_heights.append(max_h)

    # Tree panel heights
    tree_max_h = max(tp_h for _, _, tp_h in tree_panels.values()) if tree_panels else 0

    # Total image dimensions (extra ~50px for the scale legend below trees)
    SCALE_LEGEND_H = 70
    total_height = (MARGIN + sum(gp.INITIAL_Y + rh + 40 for rh in row_heights) +
                    len(rows) * ROW_GAP + tree_max_h + SCALE_LEGEND_H + MARGIN + 20)

    # Compute actual row widths accounting for multi-scaffold chroms
    # (each chrom width depends on its number of haplotypes)
    gbanding_dir = getattr(args, 'gbanding_dir', None)
    gbanding_position = getattr(args, 'gbanding_position', 'above')

    def chrom_actual_width(chrom):
        haps = chrom_draw_haps.get(chrom, [])
        if not haps:
            return 0
        base_w = len(haps) * gp.hap_width() + max(0, len(haps) - 1) * gp.HAP_GAP
        # Add overhang for short-chrom p-arm offset on either side
        extra = 0
        for hap in haps:
            scaf = chrom_map.get(chrom, {}).get(hap)
            if not scaf or scaf not in scaffold_lengths:
                continue
            scaf_len = scaffold_lengths[scaf]
            scaf_h = scaf_len * pixels_per_bp_genome
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = (cen_e - cen_s) * cen_ppb_fixed if cen_s is not None else 0
            eff_h = max(scaf_h, this_cen_h)
            if subtelo_h * 2 > eff_h:
                extra += gp.TRACK_W + gp.SUBTELO_OFFSET_GAP
        return base_w + extra

    max_row_width = 0
    for row_chroms in rows:
        drawable = [c for c in row_chroms if c in chrom_draw_haps]
        rw = sum(chrom_actual_width(c) for c in drawable)
        # Sum gaps; X/Y use the smaller HAP_GAP instead of CHROM_GAP_FIG
        for i in range(len(drawable) - 1):
            if {drawable[i], drawable[i + 1]} == {"chrX", "chrY"}:
                rw += gp.HAP_GAP
            else:
                rw += CHROM_GAP_FIG
        max_row_width = max(max_row_width, rw)
    tree_total_w = sum(tw for _, tw, _ in tree_panels.values()) + TREE_GAP * 2
    total_width = max(max_row_width, tree_total_w) + MARGIN * 2

    # Create drawing
    d = draw.Drawing(total_width, total_height, id_prefix="fig")
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # Draw chromosome rows
    print("Drawing chromosomes...")
    current_y = MARGIN

    for row_idx, row_chroms in enumerate(rows):
        draw_genome_row(
            d, row_chroms, chrom_draw_haps, chrom_map, scaffold_lengths,
            centromere_coords, all_chrom_features, colors,
            pixels_per_bp_genome, genome_height, subtelo_h, cen_ppb_fixed,
            max_cen_panel_h, args.subtelo_span, sex, background,
            x_offset=MARGIN, y_offset=current_y, chrom_gap=CHROM_GAP_FIG,
            gbanding_dir=gbanding_dir, gbanding_position=gbanding_position,
        )
        current_y += gp.INITIAL_Y + row_heights[row_idx] + 40 + ROW_GAP
        print(f"  Row {row_idx + 1}: {', '.join(row_chroms)}")

    # Draw hierarchy trees
    print("Drawing hierarchy trees...")
    tree_x = MARGIN
    tree_y = current_y + 10
    max_tree_bottom = tree_y

    for fs_name in ["chromosome", "region", "subtelomeric"]:
        if fs_name not in tree_panels:
            continue
        grp, tw, th = tree_panels[fs_name]
        container = draw.Group(transform=f"translate({tree_x},{tree_y})")
        container.append(grp)
        d.append(container)
        tree_x += tw + TREE_GAP
        max_tree_bottom = max(max_tree_bottom, tree_y + th)

    # Scale legend below the trees
    telo_ppb_fixed = subtelo_h / args.subtelo_span if args.subtelo_span > 0 else 0
    draw_scale_legend(d, MARGIN, max_tree_bottom + 20,
                      pixels_per_bp_genome, cen_ppb_fixed,
                      telo_ppb_fixed, text_color)

    save_drawing(d, args.output, total_width, total_height)


def save_drawing(d, output, total_width, total_height):
    d.save_svg(output)
    print(f"\nSaved: {output}")
    print(f"Image size: {total_width:.0f} x {total_height:.0f} px")
    out_png = Path(output).with_suffix(".png")
    from karyoplot.svg.export import svg_to_png
    if svg_to_png(output, str(out_png), dpi=300, quiet=True):
        print(f"Saved: {out_png}")


def draw_packed_layout(args, chrom_draw_haps, chrom_map, scaffold_lengths,
                       centromere_coords, all_chrom_features, colors,
                       pixels_per_bp_genome, genome_height, subtelo_h,
                       cen_ppb_fixed, max_cen_panel_h, sex, background,
                       all_chroms, tree_panels):
    """Shelf-packed chromosomes aligned by centromere, with trees in a
    contiguous column on the right (not interrupted by chromosomes)."""
    text_color = "black" if background == "white" else "white"

    CHROM_GAP_FIG = 12
    MARGIN = 30
    SHELF_GAP = 25
    HEADER_PAD = gp.INITIAL_Y  # 60
    LABEL_PAD = 45
    TREE_GAP = 20

    # -----------------------------------------------------------------------
    # Compute per-chrom element dimensions + centromere metrics
    # -----------------------------------------------------------------------
    chrom_elements = []
    drawable_chroms = [c for c in all_chroms if c in chrom_draw_haps]

    for chrom in drawable_chroms:
        hap_data = chrom_map.get(chrom, {})
        draw_haps = chrom_draw_haps[chrom]
        max_box_h = 0
        max_pre = 0    # max distance from element top to genome cen midpoint
        max_post = 0   # max distance from genome cen midpoint to element bottom
        extra_w_left = 0
        extra_w_right = 0

        for hap in draw_haps:
            scaf = hap_data.get(hap)
            if not scaf or scaf not in scaffold_lengths:
                continue
            scaf_len = scaffold_lengths[scaf]
            scaf_h = scaf_len * pixels_per_bp_genome
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            if cen_s is None:
                continue
            this_cen_h = (cen_e - cen_s) * cen_ppb_fixed
            box_h = max(scaf_h, this_cen_h)
            # Genome y_offset within box (if cen panel is taller, genome is centered)
            genome_y_off = max(0, (this_cen_h - scaf_h) / 2)
            # Centromere midpoint position in genome track
            cen_mid_bp = (cen_s + cen_e) / 2
            cen_mid_in_genome = cen_mid_bp * pixels_per_bp_genome
            # pre: from top of chrom element (header top) to cen midpoint
            pre_h = HEADER_PAD + genome_y_off + cen_mid_in_genome
            # post: from cen midpoint to bottom of chrom element
            post_h = HEADER_PAD + box_h + LABEL_PAD - pre_h
            max_pre = max(max_pre, pre_h)
            max_post = max(max_post, post_h)
            max_box_h = max(max_box_h, box_h)

            # Check if p-arm is offset (telomere overlap)
            if subtelo_h * 2 > box_h:
                is_h1 = hap.startswith("hap1")
                if sex == "male" and chrom == "chrX" and hap.startswith("hap2"):
                    is_h1 = True
                if sex == "male" and chrom == "chrY" and hap.startswith("hap1"):
                    is_h1 = False
                if is_h1:
                    extra_w_left = gp.TRACK_W + gp.SUBTELO_OFFSET_GAP
                else:
                    extra_w_right = gp.TRACK_W + gp.SUBTELO_OFFSET_GAP

        base_w = len(draw_haps) * gp.hap_width() + max(0, len(draw_haps) - 1) * gp.HAP_GAP
        elem_w = base_w + extra_w_left + extra_w_right
        elem_h = HEADER_PAD + max_box_h + LABEL_PAD

        chrom_elements.append({
            "id": chrom,
            "width": elem_w,
            "height": elem_h,
            "pre_h": max_pre,
            "post_h": max_post,
            "extra_left": extra_w_left,
        })

    # -----------------------------------------------------------------------
    # Shelf-pack chromosomes ONLY (no trees), target width ~900
    # -----------------------------------------------------------------------
    def shelf_pack(els, max_width, x_gap=CHROM_GAP_FIG):
        els_sorted = sorted(els, key=lambda e: -e["height"])
        shelves = []
        for el in els_sorted:
            placed = False
            for s in shelves:
                if s["next_x"] + el["width"] <= max_width:
                    s["items"].append(el)
                    s["next_x"] += el["width"] + x_gap
                    placed = True
                    break
            if not placed:
                shelves.append({
                    "next_x": el["width"] + x_gap,
                    "items": [el],
                })
        return shelves

    # Dynamic shelf width: at least 950, but big enough to fit ~3-4 of the
    # widest chrom elements per shelf (multi-scaffold samples have wide chroms)
    max_el_w = max(el["width"] for el in chrom_elements)
    CHROM_PACK_WIDTH = max(950, max_el_w * 4 + 50)
    shelves = shelf_pack(chrom_elements, CHROM_PACK_WIDTH)

    # Sort items within shelves by chromosome number
    for s in shelves:
        s["items"].sort(key=lambda el: gp.chromosome_sort_key(el["id"]))

    # Compute shelf dimensions using centromere alignment
    # Each chrom's cen midpoint aligns to shelf_pre_max within the shelf
    for s in shelves:
        s["pre_max"] = max(el["pre_h"] for el in s["items"])
        s["post_max"] = max(el["post_h"] for el in s["items"])
        s["height"] = s["pre_max"] + s["post_max"]
        s["width"] = sum(el["width"] for el in s["items"]) + \
                     (len(s["items"]) - 1) * CHROM_GAP_FIG

    chrom_block_w = max(s["width"] for s in shelves)
    chrom_block_h = sum(s["height"] for s in shelves) + SHELF_GAP * (len(shelves) - 1)

    # -----------------------------------------------------------------------
    # Tree block: stacked vertically on the right (contiguous)
    # -----------------------------------------------------------------------
    tree_order = ["chromosome", "region", "subtelomeric"]
    tree_items = []
    for fs_name in tree_order:
        if fs_name in tree_panels:
            grp, tw, th = tree_panels[fs_name]
            tree_items.append({"id": fs_name, "width": tw, "height": th, "group": grp})

    tree_block_w = max((t["width"] for t in tree_items), default=0)
    tree_block_h = sum(t["height"] for t in tree_items) + \
                   TREE_GAP * max(0, len(tree_items) - 1)

    # -----------------------------------------------------------------------
    # Compose: chromosomes on left, trees on right
    # -----------------------------------------------------------------------
    SCALE_LEGEND_H = 70
    total_width = chrom_block_w + 2 * MARGIN + (TREE_GAP + tree_block_w if tree_items else 0)
    total_height = max(chrom_block_h, tree_block_h) + SCALE_LEGEND_H + 2 * MARGIN

    d = draw.Drawing(total_width, total_height, id_prefix="fig")
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # Draw chromosomes
    print("Drawing chromosomes (packed + centromere-aligned)...")
    shelf_y = MARGIN
    for shelf_idx, s in enumerate(shelves):
        # Shelf centromere line (absolute y)
        shelf_cen_y = shelf_y + s["pre_max"]
        # Draw each chrom in the shelf
        cx = MARGIN
        for el in s["items"]:
            chrom = el["id"]
            # y_offset = shelf_cen_y - chrom.pre_h (so chrom cen aligns with shelf line)
            chrom_y = shelf_cen_y - el["pre_h"]
            base_x = cx + el.get("extra_left", 0)
            draw_genome_row(
                d, [chrom], chrom_draw_haps, chrom_map, scaffold_lengths,
                centromere_coords, all_chrom_features, colors,
                pixels_per_bp_genome, genome_height, subtelo_h,
                cen_ppb_fixed, max_cen_panel_h, args.subtelo_span,
                sex, background,
                x_offset=base_x, y_offset=chrom_y, chrom_gap=0,
            )
            cx += el["width"] + CHROM_GAP_FIG
        items_str = ", ".join(el["id"] for el in s["items"])
        print(f"  Shelf {shelf_idx + 1} (h={s['height']:.0f}): {items_str}")
        shelf_y += s["height"] + SHELF_GAP

    # Draw trees on right (stacked vertically)
    print("Drawing hierarchy trees (right column)...")
    tree_x = MARGIN + chrom_block_w + TREE_GAP
    tree_y = MARGIN
    for t in tree_items:
        container = draw.Group(transform=f"translate({tree_x},{tree_y})")
        container.append(t["group"])
        d.append(container)
        print(f"  {t['id']}_tree at ({tree_x}, {tree_y})")
        tree_y += t["height"] + TREE_GAP

    # Scale legend below the chromosome block
    legend_y = MARGIN + max(chrom_block_h, tree_block_h) + 20
    telo_ppb_fixed = subtelo_h / args.subtelo_span if args.subtelo_span > 0 else 0
    draw_scale_legend(d, MARGIN, legend_y,
                      pixels_per_bp_genome, cen_ppb_fixed,
                      telo_ppb_fixed, text_color)

    save_drawing(d, args.output, total_width, total_height)


if __name__ == "__main__":
    main()
