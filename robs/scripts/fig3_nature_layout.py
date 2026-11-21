#!/usr/bin/env python3
"""Re-lay the Fig. 3 SVG (Robertsonian zoom plot) to the Nature main-figure size.

Nature main figures are at most 180 mm wide and 170 mm high, with 5-7 pt Arial text
and 8 pt bold lowercase panel letters. The zoom plot is drawn with three rows of
200-unit-tall tracks (1 unit = 1 CSS px = 0.75 pt), which makes it 214 mm high.
This script shortens every track to 160 units, closes the rows up, sets the fonts,
and trims the canvas. Positions along each track are scaled linearly, so every
annotation keeps its relative position and the zoom wedges still meet the tracks.

    python3 fig3_nature_layout.py IN.svg OUT.svg

The input must be the 647 x 809 layout (track tops at y = 15, 250, 485). The legend
(x >= 470) is left where it is.
"""
import re
import sys

OLD_TRACK, NEW_TRACK = 200.0, 160.0
OLD_PITCH, NEW_PITCH = 235.0, 195.0
TOP = 15.0
N_ROWS = 3
LEGEND_X = 470.0
# Font sizes in units (0.75 pt each): old -> new. 9 units = 6.75 pt is kept.
FONT = {"16": "10.67", "13": "9.33", "12": "9.33", "10": "9.33", "9": "9"}
NEW_HEIGHT = 637  # 168.5 mm


def fy(y):
    """Map an old y coordinate to the new layout."""
    row = min(N_ROWS - 1, max(0, int((y - (TOP - 12)) // OLD_PITCH)))
    old_top = TOP + OLD_PITCH * row
    new_top = TOP + NEW_PITCH * row
    d = y - old_top
    if d <= 0:
        return new_top + d
    if d <= OLD_TRACK:
        return new_top + d * NEW_TRACK / OLD_TRACK
    return new_top + NEW_TRACK + (d - OLD_TRACK)


def num(v):
    return f"{v:.4f}".rstrip("0").rstrip(".")


def attr(line, name):
    m = re.search(rf'\b{name}="([^"]*)"', line)
    return m.group(1) if m else None


def set_attr(line, name, value):
    return re.sub(rf'\b{name}="[^"]*"', f'{name}="{value}"', line, count=1)


def do_rect(line):
    x, y, h = float(attr(line, "x")), float(attr(line, "y")), float(attr(line, "height"))
    if x >= LEGEND_X or (y < 0 and h > 700):  # legend swatch or background
        return line
    y2 = fy(y)
    line = set_attr(line, "y", num(y2))
    return set_attr(line, "height", num(fy(y + h) - y2))


def do_text(line):
    line = line.replace('font-family="sans-serif"', 'font-family="Arial"')
    size = attr(line, "font-size")
    line = set_attr(line, "font-size", FONT[size])
    if float(attr(line, "x")) >= LEGEND_X:
        return line
    y2 = num(fy(float(attr(line, "y"))))
    line = set_attr(line, "y", y2)
    return re.sub(r'rotate\((-?[\d.]+) ([\d.]+) [\d.]+\)',
                  lambda m: f"rotate({m.group(1)} {m.group(2)} {y2})", line)


def do_path(line):
    d = attr(line, "d")
    d2 = re.sub(r'([\d.]+),([\d.]+)',
                lambda m: f"{m.group(1)},{num(fy(float(m.group(2))))}", d)
    return line.replace(d, d2)


def main(src, dst):
    out = []
    for line in open(src):
        if line.lstrip().startswith("<rect"):
            line = do_rect(line)
        elif line.lstrip().startswith("<text"):
            line = do_text(line)
        elif line.lstrip().startswith("<path"):
            line = do_path(line)
        out.append(line)
    text = "".join(out)
    if 'height="809" viewBox="0 -4 647 809"' not in text:
        sys.exit("unexpected canvas: expected the 647 x 809 layout")
    text = text.replace('height="809" viewBox="0 -4 647 809"',
                        f'height="{NEW_HEIGHT}" viewBox="0 -4 647 {NEW_HEIGHT}"')
    text = text.replace('<rect x="0" y="-4" width="647" height="809" fill="white" />',
                        f'<rect x="0" y="-4" width="647" height="{NEW_HEIGHT}" fill="white" />')
    open(dst, "w").write(text)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
