#!/usr/bin/env python3
"""Assemble Supplementary Figure 5 (panels A-E) onto one canvas.

Until now this assembly existed only inside `breakpoint_density_stacked.ai`:
the committed SVG is 450x939, which is panels A-C alone, while the published
PNG carries all five. The panel letters, the relocated legend, the panel D row
labels and the leader lines were all added by hand in Illustrator and were not
reproducible. This script replaces that step.

Layout, following the published figure:

    +----------------+---------------------------+
    |  A             |  D                        |
    |  B             |                           |
    |                +---------+-----------------+
    |  C             |  legend |  E              |
    +----------------+---------+-----------------+

Panels A-C arrive as a rendered SVG (they are expensive to redraw and the
script that makes them already runs on the cluster); panels D and E are built
in-process through their `build()` functions. The legend is rebuilt here from
the feature list panels A-C emit, so it matches them exactly while sitting in
its published position rather than beside them.

The leader lines are the reason this is one script rather than three nested
SVGs: they run from panel D's row labels to the SST1 arrays in panel A, so they
cross a panel boundary and have to be drawn on the parent canvas, in a
coordinate space that knows about both.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import drawsvg as draw

sys.path.insert(0, str(Path(__file__).resolve().parent))

import chm13_kmer_lollipop as panel_e  # noqa: E402
import sst1_kmer_composition as panel_d  # noqa: E402
from breakpoint_density_stacked import draw_legend, load_colors  # noqa: E402

# ---------------------------------------------------------------------------
# Layout. Every panel keeps its own internal coordinates; these are the scale
# and position applied to each on the assembled canvas. Scales exist because the
# three panels are authored at unrelated natural sizes.
# ---------------------------------------------------------------------------
LAYOUT = {
    "margin": 10,
    "gap_col": 18,          # between the A-C column and the right-hand column
    "gap_row": 18,          # between panel D and the legend/E row
    "gap_legend_e": 18,     # between the legend and panel E
    "gap_legend_left": 12,  # legend sits this far right of the A-C column
    "scale_abc": 1.00,
    "right_width": 640,     # width of the D / (legend + E) column
    "letter_size": 15,      # the space reserved above panel E for its letter
    # Panel letters print at letter_pt when the figure is set FINAL_WIDTH_MM wide
    # (Nature Extended Data: 8 pt bold lowercase at 183 mm); the size in canvas
    # units follows from the derived canvas width.
    "letter_pt": 8,
    # scale_d, scale_e and scale_legend are DERIVED, not set. The three panels
    # are authored at unrelated natural sizes (A-C 338x969, D 395x417,
    # E 353x565), so fixed scales leave the right-hand column a different height
    # from the A-C column and the figure ends up with dead space down one side.
    # Panel D is scaled to fill the column width, then the legend/E row is
    # scaled to consume exactly the height A-C leaves over.
}

# Panel A is GM03417; the published figure draws leader lines only for that
# sample, from its three SST1 arrays out to panel D's three row labels.
FINAL_WIDTH_MM = 183

LEADER_SAMPLE = "GM03417"
LEADER_COLOR = "#000000"
LEADER_WIDTH = 0.7


def embed_svg(path, x, y, scale):
    """Nest an external SVG at (x, y), scaled. Returns (raw, width, height)."""
    txt = Path(path).read_text()
    m = re.search(r'<svg[^>]*viewBox="([^"]+)"[^>]*>(.*)</svg>', txt, re.S)
    if not m:
        raise SystemExit(f"could not parse an SVG viewBox out of {path}")
    vb, inner = m.group(1), m.group(2)
    _, _, vw, vh = (float(v) for v in vb.replace(",", " ").split())
    w, h = vw * scale, vh * scale
    raw = draw.Raw(
        f'<svg x="{x}" y="{y}" width="{w}" height="{h}" viewBox="{vb}">'
        f"{inner}</svg>")
    return raw, w, h


def panel_letter(d, letter, x, y, size):
    d.append(draw.Text(letter, font_size=size, x=x, y=y, fill="#000000",
                       font_family="sans-serif", font_weight="bold",
                       text_anchor="start", dominant_baseline="hanging"))


def sst1_array_y(anchors, sample, arrays, scale, y_off):
    """Assembled-canvas y for the centre of each of `arrays`, in panel A.

    Anchored on panel D's own array intervals rather than on the drawn SST1
    marks. The marks fragment into dozens of sub-pixel runs (dozens per sample),
    so clustering them would be guesswork; panel D's coordinates are the
    definition of what each row means, which is what the line should point at.
    """
    row = next((s for s in anchors["samples"] if s["sample"] == sample), None)
    if row is None:
        raise SystemExit(f"no anchor record for {sample}")
    span = row["focus_end"] - row["focus_start"]
    ppbp = row["panel_height"] / span
    out = []
    for _label, start, end in arrays:
        mid = (start + end) / 2
        y_local = row["panel_y"] + (mid - row["focus_start"]) * ppbp
        out.append(y_off + y_local * scale)
    x_local = row["sst1_x"] + row["sst1_width"]
    return out, x_local


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--abc-svg", required=True,
                    help="Rendered panels A-C (use --no-legend when making it).")
    ap.add_argument("--abc-anchors", required=True,
                    help="Anchor JSON from --emit-anchors.")
    ap.add_argument("--bed-dir", required=True, help="Panel D's tabix BEDs.")
    ap.add_argument("--colors-dir", required=True)
    ap.add_argument("--region-bed", required=True, help="Panel E's v1 CHM13 BED.")
    ap.add_argument("--chroms", default="chr7,chr9,chr17,chr20",
                    help="Panel E chromosomes. The published figure omits chr8, "
                         "which carries a single k-mer.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    L = LAYOUT
    anchors = json.loads(Path(args.abc_anchors).read_text())

    # ---- panels ---------------------------------------------------------
    print("Building panel D...")
    panel_d.load_feature_colors(args.colors_dir)
    d_data, d_present = panel_d.collect_data(args.bed_dir, verbose=False)
    g_d, d_w, d_h, d_anchors = panel_d.build(d_data, d_present, draw_title=False)

    # Panel E is built AFTER the row-2 height is known, so it can be authored at
    # that height rather than scaled down into it -- a transform would shrink its
    # chromosome labels and lollipop counts with it, and they are already small.
    print("Building panel E...")
    g_e = e_w = e_h = None  # built below, once row2_h is known

    # ---- placement ------------------------------------------------------
    m = L["margin"]
    abc_x, abc_y = m, m
    raw_abc, abc_w, abc_h = embed_svg(args.abc_svg, abc_x, abc_y, L["scale_abc"])

    right_x = abc_x + abc_w + L["gap_col"]

    # ---- the alignment grid --------------------------------------------
    # Shared edges, so the panel rectangles line up across the whole figure:
    #
    #   y = m                top of A          = top of D
    #   y = box_B_bot        bottom of B       = bottom of D (its "% of k-mers")
    #   y = box_C_top        top of C          = top of legend = top of E
    #   y = box_C_bot        bottom of C       = bottom of legend
    #   y = abc_bottom       bottom of A-C incl. track labels = bottom of E
    #
    #   x = right_x          left of D         = left of legend
    #   x = right_x+right_w  right of D        = right of E
    #
    # Panel D is sized by WIDTH here (it fills the column) and its bottom is
    # reported against box_B_bot; the panel's own bar pitch is tuned so the two
    # land together. Panel E runs deeper than the legend on purpose -- that is
    # what makes its chromosomes long, and it fills what would be dead space
    # beside the A-C track labels.
    boxes = [(abc_y + r["panel_y"] * L["scale_abc"],
              abc_y + (r["panel_y"] + r["panel_height"]) * L["scale_abc"])
             for r in anchors["samples"]]
    box_A_top = boxes[0][0]
    box_B_bot = boxes[1][1]
    box_C_top, box_C_bot = boxes[2]
    abc_bottom = abc_y + abc_h

    ks_colors = load_colors(args.colors_dir, "KS_human_CHM13",
                            ["chromosome", "acrocentric", "region", "repeat", "gene"],
                            background="white")
    displayed = {k: set(v) for k, v in anchors["displayed"].items()}
    g_legend = draw.Group()
    # draw_legend sizes its box from the MEASURED text, so take the width it
    # reports rather than assuming one -- an assumed 236 px is what let the two
    # longest Density Analysis labels cross the border under rsvg-convert.
    legend_h, legend_box_w = draw_legend(g_legend, 0, 0, ks_colors,
                                         ["chromosome", "acrocentric", "region", "repeat", "gene"],
                                         displayed, "#000000")
    legend_w = legend_box_w + 10   # published breathing room between the border and E
    scale_legend = (box_C_bot - box_C_top) / legend_h
    legend_sw = legend_w * scale_legend

    # Panel E is authored at its slot height so its text keeps full size.
    # 38 px of that height is fixed chrome (chromosome labels above, scale-bar
    # caption below).
    # Lift E so the TOP of its chromosome labels lines up with the top of the
    # letter E, as the published figure has it. The labels sit ~5 px below E's
    # own origin (margin_top 22, baseline at -8, cap height ~9), and the letter's
    # top is ls+2 above box_C_top.
    e_lift = L["letter_size"] + 2 + 5
    e_y_top = box_C_top - e_lift
    e_slot_h = abc_bottom - e_y_top
    g_e, e_w, e_h, _ = panel_e.build(args.region_bed,
                                     chroms=args.chroms.split(","),
                                     draw_title=False,
                                     track_height=max(80.0, e_slot_h - 38))
    scale_e = 1.0
    e_sw, e_sh = e_w, e_h

    # The legend sits LEFT of panel D, as in the published figure, rather than
    # sharing D's left edge -- so the right-hand column is measured from where
    # the legend actually starts.
    legend_x = abc_x + abc_w + L["gap_legend_left"]
    # The column width is set so that panel D, scaled to fill it, ends on panel
    # B's bottom edge. The legend's width depends on the font's measured text
    # widths, so the gap between the legend and E takes up the difference
    # (it is gap_legend_e when the font is the one the layout was tuned with).
    right_w = (box_B_bot - m) / d_h * d_w
    gap_legend_e = right_w - (legend_x - right_x) - legend_sw - e_sw
    if gap_legend_e < L["gap_legend_e"]:
        right_w += L["gap_legend_e"] - gap_legend_e
        gap_legend_e = L["gap_legend_e"]
    print(f"  legend-E gap {gap_legend_e:.1f} px")
    scale_d = right_w / d_w
    d_sw, d_sh = right_w, d_h * scale_d
    d_x, d_y = right_x, m

    d_bottom = d_y + d_sh
    print(f"  grid: D bottom {d_bottom:.0f} vs panel B bottom {box_B_bot:.0f} "
          f"(delta {d_bottom - box_B_bot:+.0f} px)")

    # draw_legend() draws its border at (x-5, y-5), so placing the group at
    # box_C_top put the visible RECTANGLE 5*scale above the line. Shifting down
    # by exactly that puts the border's top and bottom on panel C's.
    legend_y = box_C_top + 5 * scale_legend
    e_x, e_y = right_x + right_w - e_sw, e_y_top
    row2_h = e_slot_h

    total_w = right_x + right_w + m
    total_h = (abc_y + abc_h) + m

    # ---- draw -----------------------------------------------------------
    d = draw.Drawing(total_w, total_h)
    d.append(draw.Rectangle(0, 0, total_w, total_h, fill="white"))
    d.append(raw_abc)

    d.append(draw.Group(transform=f"translate({d_x},{d_y}) scale({scale_d})",
                        children=[g_d]))
    d.append(draw.Group(transform=f"translate({e_x},{e_y}) scale({scale_e})",
                        children=[g_e]))
    d.append(draw.Group(
        transform=f"translate({legend_x},{legend_y}) scale({scale_legend})",
        children=[g_legend]))

    # ---- leader lines: panel A's SST1 arrays -> panel D's row labels -----
    ys, sst1_x_local = sst1_array_y(anchors, LEADER_SAMPLE,
                                    panel_d.SAMPLES[0]["arrays"],
                                    L["scale_abc"], abc_y)
    x_from = abc_x + sst1_x_local * L["scale_abc"]
    for i, (lx, ly) in enumerate(d_anchors["row_label_anchors"]):
        if i >= len(ys):
            break
        # Stop clear of the rotated label rather than crossing its glyphs.
        x_to = d_x + (lx - 9) * scale_d
        y_to = d_y + ly * scale_d
        d.append(draw.Line(x_from, ys[i], x_to, y_to,
                           stroke=LEADER_COLOR, stroke_width=LEADER_WIDTH,
                           stroke_opacity=0.85))

    # ---- panel letters --------------------------------------------------
    ls = round(L["letter_pt"] * total_w / (FINAL_WIDTH_MM / 25.4 * 72), 2)
    for letter, (top, _bot) in zip("abc", boxes):
        panel_letter(d, letter, abc_x + 4, top - ls - 2, ls)
    # D sits level with A, E level with C, as requested.
    panel_letter(d, "d", d_x - 2, box_A_top - ls - 2, ls)
    panel_letter(d, "e", e_x - 15, box_C_top - ls - 2, ls)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    d.save_svg(str(out))
    print(f"Saved {out}  ({total_w:.0f} x {total_h:.0f})")
    print(f"  aspect w/h {total_w / total_h:.3f}  (manuscript 0.767)")
    print(f"  A-C {abc_w:.0f}x{abc_h:.0f}  D {d_sw:.0f}x{d_sh:.0f} (x{scale_d:.2f})  "
          f"E {e_sw:.0f}x{e_sh:.0f}  legend x{scale_legend:.2f}")

    import subprocess
    png = out.with_suffix(".png")
    try:
        # FINAL_WIDTH_MM at 300 ppi, on white. "-d 300" alone sets only the DPI
        # tag and gives a 1:1 raster (~870 px).
        width_px = str(round(FINAL_WIDTH_MM / 25.4 * 300))
        subprocess.run(["rsvg-convert", "-w", width_px, "-b", "white",
                        str(out), "-o", str(png)], check=True, capture_output=True)
        print(f"Saved {png}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Note: rsvg-convert not available; PNG not generated")


if __name__ == "__main__":
    main()
