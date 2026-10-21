#!/usr/bin/env python3
"""Add 4 units of headroom to the rendered Fig. 3 SVG, as in the submitted figure.

The conductor render puts the top scale labels slightly above the canvas. This edit
extends the canvas upward by PAD units: the viewBox starts at -PAD, the height grows
by PAD, and the white background rect starts at y = -PAD. Nothing else changes.

    python3 fig3_headroom.py robertsonian_integrated.svg out.svg
"""
import re
import sys

PAD = 4


def main(src, dst):
    s = open(src).read()
    m = re.search(r'width="(\d+)" height="(\d+)" viewBox="0 0 (\d+) (\d+)"', s)
    if not m:
        sys.exit("fig3_headroom: expected an SVG header with viewBox starting at 0 0")
    w, h = m.group(1), int(m.group(2))
    s = s.replace(m.group(0), f'width="{w}" height="{h + PAD}" viewBox="0 -{PAD} {w} {h + PAD}"', 1)
    old = f'<rect x="0" y="0" width="{w}" height="{h}" fill="white" />'
    if old not in s:
        sys.exit("fig3_headroom: background rect not found")
    s = s.replace(old, f'<rect x="0" y="-{PAD}" width="{w}" height="{h + PAD}" fill="white" />', 1)
    open(dst, "w").write(s)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
