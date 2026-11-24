#!/usr/bin/env python3
"""Draw featureset hierarchy trees for Fig. 1 at final size, one SVG per tree.

The trees are those of Extended Data Figs. 1b, 2b and 3b: plot_featureset_tree_bars.draw_panel on
the v2 hierarchy (KS_human_CHM13_v2.priority.txt, the k-mer counts of the database info files, the
colour files), with its geometry, and with text sized in printed points:

    python draw_fig1_trees.py --featuresets region subtelomeric --label-pt 5.5 --out-dir trees/

writes trees/tree_region.svg and trees/tree_subtelomeric.svg in pt units (1 SVG unit = 1 pt) and
prints each size in mm. The geometry is scaled so the labels print at --label-pt; the root label
(the header), the k-mer counts and the axis ticks are set at their own sizes. --max-bar-w widens the
k-mer bars so the axis tick labels (1, 1.0K, 1.0M, 1.0B) do not overlap at those sizes. Text widths
are Arial advances. KS_RESOURCES must point to vendored/KaryoScope-archive-resources/resources.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import drawsvg as draw
from PIL import ImageFont

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "vendored" / "KaryoScope-ISCN"))
import plot_featureset_tree_bars as tp  # noqa: E402

MM = 72 / 25.4


def arial(bold):
    path = subprocess.run(["fc-match", "-f", "%{file}", "Arial:bold" if bold else "Arial"],
                          capture_output=True, text=True, check=True).stdout.strip()
    return ImageFont.truetype(path, 1000)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--featuresets", nargs="+", default=["region"])
    ap.add_argument("--label-pt", type=float, default=5.5, help="node labels (pt)")
    ap.add_argument("--root-pt", type=float, default=None, help="root label, the tree's header (default: label + 1)")
    ap.add_argument("--count-pt", type=float, default=5.0, help="k-mer counts (pt)")
    ap.add_argument("--tick-pt", type=float, default=5.0, help="axis tick labels (pt)")
    ap.add_argument("--max-bar-w", type=float, default=None,
                    help="longest k-mer bar, in the generator's units (60 there); default: the shortest "
                         "at which the axis tick labels are --tick-gap-pt apart")
    ap.add_argument("--tick-gap-pt", type=float, default=1.0)
    ap.add_argument("--width-mm", type=float, default=None,
                    help="make each tree this wide by setting --max-bar-w (the longest k-mer bar)")
    ap.add_argument("--no-header", action="store_true",
                    help="leave out the bold root label (the feature set's name) above the tree")
    ap.add_argument("--axis-title", default=None,
                    help='title of the k-mer axis, set on the tick-label row left of the first tick; '
                         'text between asterisks is italic, e.g. "Number of *k*-mers"')
    ap.add_argument("--axis-title-pt", type=float, default=None, help="axis title (pt; default: --tick-pt)")
    ap.add_argument("--no-multigroup", action="store_true", help='leave out " (multigroup)" after inner-node counts')
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    root_pt = a.root_pt or a.label_pt + 1.0
    f = a.label_pt / tp.FONT_SIZE                       # pt per generator unit
    fonts = {False: arial(False), True: arial(True)}

    def text_width(text, size, bold):
        return fonts[bool(bold)].getlength(text) / 1000 * size

    trees = tp.load_all_featuresets_v2()
    for root in trees.values():
        tp._apply_color_overrides(root)
    info = {}
    for fs in a.featuresets:
        min_bar = tp.min_bar_w_for_ticks(trees[fs], a.tick_pt / f, text_width, a.tick_gap_pt / f)
        bar = a.max_bar_w or min_bar
        if bar < min_bar - 1e-9:
            sys.exit(f"{fs}: --max-bar-w {bar} lets the axis tick labels touch (needs {min_bar:.1f})")
        style = dict(font_size=tp.FONT_SIZE, root_font=root_pt / f, count_font=a.count_pt / f,
                     tick_font=a.tick_pt / f, max_bar_w=bar, font_family="Arial",
                     text_width=text_width, multigroup=not a.no_multigroup, min_w=0)
        if a.no_header:
            style["root_label"] = False
        if a.axis_title:
            parts = [(t, i % 2 == 1) for i, t in enumerate(a.axis_title.split("*")) if t]
            style["axis_title"] = parts
            style["axis_title_font"] = (a.axis_title_pt or a.tick_pt) / f
        g, w, h = tp.draw_panel(trees[fs], fs, style)
        if a.width_mm:                     # width grows with the bar length: bisect for it
            lo, hi = min_bar, 4000.0
            for _ in range(60):
                style["max_bar_w"] = (lo + hi) / 2
                g, w, h = tp.draw_panel(trees[fs], fs, style)
                lo, hi = (style["max_bar_w"], hi) if w * f < a.width_mm * MM else (lo, style["max_bar_w"])
            print(f"  {fs}: longest bar {style['max_bar_w'] * f / MM:.1f} mm for a {a.width_mm} mm tree")
        # the axis tick labels sit above the root row: room for them above the group
        top = max(0.0, a.tick_pt / f + 5 - (4 + tp.ROW_H / 2) + 1)
        if a.axis_title:                   # a title taller than the ticks needs its own room
            top = max(top, 0.75 * (a.axis_title_pt or a.tick_pt) / f - (4 + tp.ROW_H / 2 - 5) + 0.5)
        # right edge: draw_panel's estimate plus nothing; measured widths are exact
        W, H = w * f, (h + top) * f
        d = draw.Drawing(W, H)
        d.set_render_size(f"{W:.4f}pt", f"{H:.4f}pt")    # 1 unit = 1 pt (unitless would be px)
        d.append(draw.Group([g], transform=f"scale({f}) translate(0,{top})"))
        p = out / f"tree_{fs}.svg"
        d.save_svg(str(p))
        info[fs] = dict(width_mm=round(W / MM, 2), height_mm=round(H / MM, 2), scale=f, top_units=top)
        print(f"{p}: {W / MM:.1f} x {H / MM:.1f} mm (longest bar {style['max_bar_w'] * f / MM:.1f} mm; labels {a.label_pt} pt, "
              f"header {'none' if a.no_header else f'{root_pt} pt'}, counts {a.count_pt} pt, ticks {a.tick_pt} pt"
              + (f", axis title {a.axis_title!r} {a.axis_title_pt or a.tick_pt} pt" if a.axis_title else "") + ")")
    (out / "trees.json").write_text(json.dumps(info, indent=1))


if __name__ == "__main__":
    sys.exit(main())
