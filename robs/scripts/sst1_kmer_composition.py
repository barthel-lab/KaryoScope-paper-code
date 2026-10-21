#!/usr/bin/env python3
"""Bar plot of chromosome-specific k-mer composition at SST1 arrays (drawsvg)."""

import math
import subprocess
import numpy as np
import drawsvg as draw
from pathlib import Path

# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Samples and SST1 array coordinates
# ---------------------------------------------------------------------------
def _load_sst1_samples():
    """Read the SST1 array table from the `sst1_arrays` section of the ROBs config.

    The config path is $ROBS_CONFIG, or conf/robs.yaml under the project root.
    Each entry: {sample, asm, scaffold, iscn, arrays: [[label, start, end], ...]}.
    """
    import os
    import yaml
    path = Path(os.environ.get("ROBS_CONFIG", PROJECT_ROOT / "conf" / "robs.yaml"))
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    out = []
    for e in cfg["sst1_arrays"]:
        out.append({
            "sample": e["sample"], "asm": e["asm"], "scaffold": e["scaffold"],
            "iscn": e["iscn"],
            "arrays": [(a[0], int(a[1]), int(a[2])) for a in e["arrays"]],
        })
    return out


SAMPLES = _load_sst1_samples()

# Feature display order (chromosome-specific only, no multigroup/novel/other)
FEATURE_ORDER = [
    "chr13_specific",
    "chr14_specific",
    "chr21_specific",
    "chr9_specific",
    "chr20_specific",
    "chr17_specific",
    "chr7_specific",
]

FEATURE_LABELS = {
    "chr13_specific": "chr13",
    "chr14_specific": "chr14",
    "chr21_specific": "chr21",
    "chr9_specific": "chr9",
    "chr20_specific": "chr20",
    "chr17_specific": "chr17",
    "chr7_specific": "chr7",
}

# ---------------------------------------------------------------------------
# Load colors from KaryoScope chromosome color file
# ---------------------------------------------------------------------------
FEATURE_COLORS = {}

# The published row labels. The code's internal names describe POSITION relative
# to the main array; these describe what the blocks are. Note the flanking rows
# deliberately share one label -- they are the same feature type on either side,
# told apart by their row, exactly as in the manuscript figure.
ARRAY_DISPLAY_LABELS = {
    "Upstream partial": "2S/2L partial monomer",
    "Main array": "SST1 sub family 1 array",
    "Downstream partial": "2S/2L partial monomer",
}


def load_feature_colors(colors_dir):
    """Load chromosome colours, failing loudly if they are not there.

    The old module-level version pointed at a macOS `~/Documents/...` path and
    silently left FEATURE_COLORS empty when it did not exist, so every bar fell
    back to grey #545454 and the figure rendered wrong at exit 0.
    """
    path = Path(colors_dir) / "KS_human_CHM13.chromosome.colors.txt"
    if not path.exists():
        raise SystemExit(
            f"chromosome colours not found: {path}\n"
            "Pass --colors-dir pointing at the directory holding "
            "KS_human_CHM13.chromosome.colors.txt. Without it every bar would "
            "render grey and the script would still exit 0.")
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] != "feature":
                FEATURE_COLORS[parts[0]] = parts[1]
    if not FEATURE_COLORS:
        raise SystemExit(f"no colours parsed from {path}")
    return FEATURE_COLORS


def _get_color(feat):
    if feat in FEATURE_COLORS:
        return FEATURE_COLORS[feat]
    return "#545454"


# ---------------------------------------------------------------------------
# Query chromosome featureset via tabix
# ---------------------------------------------------------------------------
def query_composition(asm, scaffold, start, end, bed_dir=None):
    bed_dir = Path(bed_dir) if bed_dir else (PROJECT_ROOT / "data" / "robs")
    bed_path = (bed_dir /
                f"{asm}.KS_human_CHM13.chromosome.presmoothed.scaffolded.features.bed.gz")
    if not bed_path.exists():
        raise SystemExit(
            f"chromosome BED not found: {bed_path}\n"
            "This must be BGZF-compressed and tabix-indexed (the staged "
            "remap-bed output is plain gzip, which tabix cannot read).")
    if not Path(str(bed_path) + ".tbi").exists():
        raise SystemExit(f"missing tabix index: {bed_path}.tbi  (run: tabix -p bed <file>)")
    region = f"{scaffold}:{start}-{end}"
    result = subprocess.run(
        ["tabix", str(bed_path), region],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"tabix failed on {region}: {result.stderr.strip()}")
    if not result.stdout.strip():
        # An empty region yields empty bars rather than an error, so raise.
        raise SystemExit(
            f"tabix returned NO records for {region} in {bed_path.name}. "
            "Check the scaffold name matches the scaffolded frame.")
    bp_by_feat = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        feat = parts[3]
        s, e = int(parts[1]), int(parts[2])
        os_ = max(s, start)
        oe = min(e, end)
        if oe > os_:
            bp_by_feat[feat] = bp_by_feat.get(feat, 0) + (oe - os_)

    total = sum(bp_by_feat.values())
    if total == 0:
        return {}, 0

    known = set(FEATURE_ORDER)
    pct = {}
    for feat, bp in bp_by_feat.items():
        if feat in known:
            pct[feat] = bp / total * 100
    return pct, total


# ---------------------------------------------------------------------------
# Build data
# ---------------------------------------------------------------------------
def collect_data(bed_dir=None, verbose=True):
    """Query composition for every sample x array. Returns (all_data, present)."""
    if verbose:
        print("Querying chromosome featureset composition at SST1 arrays...")
    all_data = []
    for s in SAMPLES:
        sample_data = []
        for label, start, end in s["arrays"]:
            pct, total = query_composition(s["asm"], s["scaffold"], start, end, bed_dir)
            sample_data.append({"label": label, "size_kb": (end - start) / 1000,
                                "pct": pct})
            if verbose:
                print(f"  {s['sample']} {label}: {total} bp, {len(pct)} features")
        all_data.append(sample_data)
    present = [f for f in FEATURE_ORDER
               if any(arr["pct"].get(f, 0) > 0 for sd in all_data for arr in sd)]
    if not present:
        raise SystemExit("no chromosome-specific features present in any array")
    return all_data, present


# ---------------------------------------------------------------------------
# Build the panel
# ---------------------------------------------------------------------------
def build(all_data, present_feats, draw_title=False):
    """Return (group, width, height, anchors) for panel D.

    The group is in its own coordinate space starting at (0, 0); the assembler
    translates it. `anchors` carries what the assembler cannot recompute:

      row_label_anchors  one (x, y) per array row, at the LEFT edge of that
                         row's label column -- where the leader lines to panel
                         A's SST1 arrays terminate
      row_labels         the display label for each row
      rows_y             (top, bottom) of each row band

    `draw_title` is off by default: the standalone render wants a title, the
    assembled figure gets its panel letter instead.
    """
    n_feats = len(present_feats)
    font_family = "sans-serif"
    font_size_title = 10
    font_size_label = 10
    font_size_axis = 10    # the 0.01%/0.1%/1%/10% ticks
    font_size_tick = 8
    text_color = "#000000"
    axis_color = "#545454"

    # Widened from 78 to keep the panel's aspect at ~0.79 after the taller
    # caption margin; that is what lands D's bottom on panel B's bottom line
    # in the assembled figure.
    bar_axis_len = 88
    feat_label_w = 40
    panel_w = feat_label_w + bar_axis_len
    # 12 made this panel landscape (410x407); the published panel D is
    # portrait (aspect 0.67) with visibly thicker bars. 18 gives 410x518,
    # which lets it fill the right column's width AND end level with panel B.
    bar_group_h_total = n_feats * 18
    panel_h = bar_group_h_total + 10
    bar_h = 13
    bar_y_pad = 5

    n_cols = len(SAMPLES)
    n_rows = 3
    col_gap = 4         # samples sit close; the boxes separate them
    row_gap = 14
    margin_left = 20    # column for the rotated row labels
    margin_top = 55 if draw_title else 40
    margin_bottom = 56   # rotated tick labels (~34 px) + the axis caption
    margin_right = 10

    total_w = margin_left + n_cols * panel_w + (n_cols - 1) * col_gap + margin_right
    total_h = margin_top + n_rows * panel_h + (n_rows - 1) * row_gap + margin_bottom

    x_min_pct, x_max_pct = 0.008, 30.0
    log_min, log_max = math.log10(x_min_pct), math.log10(x_max_pct)
    log_range = log_max - log_min

    def pct_to_x(pct_val, bar_origin_x):
        if pct_val <= 0:
            return bar_origin_x
        frac = (math.log10(max(pct_val, x_min_pct)) - log_min) / log_range
        return bar_origin_x + bar_axis_len * frac

    g = draw.Group()

    if draw_title:
        g.append(draw.Text(
            "Chromosome-specific k-mer composition at SST1 arrays",
            font_size=font_size_title, x=total_w / 2, y=14,
            fill=text_color, font_family=font_family,
            text_anchor="middle", font_weight="bold"))

    rows_y = []
    for si, (sample_info, sample_data) in enumerate(zip(SAMPLES, all_data)):
        for ai, arr in enumerate(sample_data):
            px = margin_left + si * (panel_w + col_gap)
            py = margin_top + ai * (panel_h + row_gap)
            bar_origin_x = px + feat_label_w
            if si == 0:
                rows_y.append((py, py + panel_h))

            g.append(draw.Rectangle(bar_origin_x, py, bar_axis_len, panel_h,
                                    fill="white", stroke=axis_color, stroke_width=0.5))

            if ai == 0:
                # Centred on the plot box. px + panel_w/2 centred on the box
                # PLUS its left label gutter, which put every header ~20 px left
                # of the column it names.
                g.append(draw.Text(
                    sample_info["sample"], font_size=font_size_tick + 2,
                    x=bar_origin_x + bar_axis_len / 2, y=py - 12,
                    fill=text_color, font_family=font_family,
                    text_anchor="middle", font_weight="bold"))

            # Array size, inside the box as in the published figure. Top row
            # puts it at the top, the lower two at the bottom -- following the
            # published figure, and it is also where the bars are not.
            if ai == 0:
                sy, sbase = py + 3, "hanging"
            else:
                sy, sbase = py + panel_h - 3, "auto"
            g.append(draw.Text(
                f"{arr['size_kb']:.1f} kb", font_size=font_size_tick,
                x=bar_origin_x + bar_axis_len - 3, y=sy,
                fill=text_color, font_family=font_family,
                text_anchor="end", dominant_baseline=sbase))

            for tv in (0.01, 0.1, 1, 10):
                tx = pct_to_x(tv, bar_origin_x)
                if not (bar_origin_x <= tx <= bar_origin_x + bar_axis_len):
                    continue
                g.append(draw.Line(tx, py, tx, py + panel_h, stroke=axis_color,
                                   stroke_width=0.3, stroke_opacity=0.3))
                if ai == n_rows - 1:
                    lbl = f"{tv:.0f}%" if tv >= 1 else (f"{tv:.1f}%" if tv >= 0.1
                                                       else f"{tv:.2f}%")
                    # Rotated: at this column width the horizontal labels
                    # ("0.01%", "0.1%", "1%", "10%") ran into one another.
                    ty = py + panel_h + 4
                    g.append(draw.Text(lbl, font_size=font_size_axis,
                                       x=tx, y=ty, fill=axis_color,
                                       font_family=font_family, text_anchor="end",
                                       dominant_baseline="middle",
                                       transform=f"rotate(-90 {tx} {ty})"))

            for fi, feat in enumerate(present_feats):
                by = py + bar_y_pad + fi * (bar_group_h_total / n_feats)
                bar_center_y = by + (bar_group_h_total / n_feats) / 2
                if si == 0:
                    g.append(draw.Text(
                        FEATURE_LABELS[feat], font_size=font_size_tick,
                        x=bar_origin_x - 3, y=bar_center_y, fill=_get_color(feat),
                        font_family=font_family, text_anchor="end",
                        dominant_baseline="middle"))
                val = arr["pct"].get(feat, 0)
                if val <= 0:
                    continue
                bw = max(pct_to_x(val, bar_origin_x) - bar_origin_x, 0.5)
                g.append(draw.Rectangle(bar_origin_x, by + 1, bw, bar_h,
                                        fill=_get_color(feat)))

        if si == 1:
            g.append(draw.Text(
                "% of k-mers", font_size=font_size_label,
                x=margin_left + si * (panel_w + col_gap) + feat_label_w
                  + bar_axis_len / 2,
                y=total_h - 6, fill=text_color, font_family=font_family,
                text_anchor="middle"))

    # Rotated row labels down the left edge, one per array row.
    # Row-label font is set from the row height: at font_size 10 the longest
    # label ("SST1 sub family 1 array") is ~120 px against a ~94 px row, so it
    # would overflow into the neighbouring rows and they would collide.
    row_h = rows_y[0][1] - rows_y[0][0]
    longest = max(len(ARRAY_DISPLAY_LABELS.get(a[0], a[0]))
                  for a in SAMPLES[0]["arrays"])
    font_size_row = min(font_size_label, max(6.0, (row_h - 4) / (longest * 0.55)))
    row_label_x = margin_left - 12
    row_labels, row_label_anchors = [], []
    for ai, (top, bot) in enumerate(rows_y):
        internal = SAMPLES[0]["arrays"][ai][0]
        label = ARRAY_DISPLAY_LABELS.get(internal, internal)
        row_labels.append(label)
        cy = (top + bot) / 2
        g.append(draw.Text(label, font_size=font_size_row,
                           x=row_label_x, y=cy, fill=text_color,
                           font_family=font_family, text_anchor="middle",
                           dominant_baseline="auto",
                           transform=f"rotate(-90 {row_label_x} {cy})"))
        row_label_anchors.append((row_label_x, cy))

    anchors = {"row_labels": row_labels,
               "row_label_anchors": row_label_anchors,
               "rows_y": rows_y}
    return g, total_w, total_h, anchors


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed-dir", default=None,
                    help="Directory holding the BGZF+tabix chromosome BEDs.")
    ap.add_argument("--colors-dir", required=True,
                    help="Directory with KS_human_CHM13.chromosome.colors.txt")
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    load_feature_colors(args.colors_dir)
    all_data, present = collect_data(args.bed_dir)
    g, w, h, anchors = build(all_data, present, draw_title=True)

    d = draw.Drawing(w, h)
    d.append(draw.Rectangle(0, 0, w, h, fill="white"))
    d.append(g)

    out_dir = Path(args.output_dir) if args.output_dir else (
        PROJECT_ROOT / "results" / "breakpoint_density")
    out_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / "sst1_kmer_composition.svg"
    d.save_svg(str(svg_path))
    print(f"Saved {svg_path}  ({w:.0f} x {h:.0f})")
    print("Row labels:", anchors["row_labels"])

    png_path = out_dir / "sst1_kmer_composition.png"
    try:
        subprocess.run(["rsvg-convert", "-d", "300", "-p", "300",
                        str(svg_path), "-o", str(png_path)],
                       check=True, capture_output=True)
        print(f"Saved {png_path}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Note: rsvg-convert not available; PNG not generated")


if __name__ == "__main__":
    main()
