#!/usr/bin/env python3
"""The RGB colours of the sources of each Illustrator panel, for ai_cmyk_to_rgb.py.

    python source_palettes.py --svg karyo_centromere.region=HG002...centromere.region.KaryoScope.svg \\
        --tree region --tree repeat --tree gene --tree acrocentric --tree subtelomeric --out palettes.json

writes {panel name: {hex colour: count}}: for --svg NAME=FILE every fill and stroke colour of the SVG,
for --tree FS ("tree_FS") every colour of the featureset tree as plot_featureset_tree_bars.draw_panel
draws it (v2 hierarchy, with its colour overrides). KS_RESOURCES must point to
vendored/KaryoScope-archive-resources/resources.
"""
import argparse
import collections
import json
import re
import sys
from pathlib import Path

import drawsvg as draw

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "vendored" / "KaryoScope-ISCN"))
import plot_featureset_tree_bars as tp  # noqa: E402


def colours(svg_text):
    return dict(collections.Counter(h.lower() for h in re.findall(r'(?:fill|stroke)\s*[=:]\s*"?(#[0-9A-Fa-f]{6})', svg_text)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--svg", action="append", default=[], metavar="NAME=FILE")
    ap.add_argument("--tree", action="append", default=[], metavar="FEATURESET")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pal = {}
    for spec in a.svg:
        name, path = spec.split("=", 1)
        pal[name] = colours(Path(path).read_text())
    if a.tree:
        trees = tp.load_all_featuresets_v2()
        for root in trees.values():
            tp._apply_color_overrides(root)
        for fs in a.tree:
            g, w, h = tp.draw_panel(trees[fs], fs)
            d = draw.Drawing(w, h)
            d.append(g)
            pal["tree_" + fs] = colours(d.as_svg())
    Path(a.out).write_text(json.dumps(pal, indent=1))
    for k, v in pal.items():
        print(f"{k}: {len(v)} colours")


if __name__ == "__main__":
    sys.exit(main())
