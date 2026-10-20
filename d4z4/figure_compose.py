#!/usr/bin/env python3
"""
Compose a multi-panel figure page from panel SVGs.

A page is a list of SVG fragments in page points (1/72 inch, y down): panel SVGs placed with
nest_svg(), plus text and lines. Panels are inlined (their <svg> root becomes a nested <svg> with x/y/width/height and the
original viewBox), so text stays text and every element stays editable in Illustrator.
"""

import re
from pathlib import Path
from xml.sax.saxutils import escape

def _strip_prolog(svg_text):
    svg_text = re.sub(r'<\?xml[^>]*\?>', '', svg_text)
    svg_text = re.sub(r'<!DOCTYPE[^>]*>', '', svg_text)
    return svg_text.strip()


def _prefix_ids(svg_text, prefix):
    """Make element ids unique per panel (matplotlib and drawsvg both reuse id names)."""
    svg_text = re.sub(r'\bid="([^"]+)"', lambda m: f'id="{prefix}{m.group(1)}"', svg_text)
    svg_text = re.sub(r'url\(#([^)]+)\)', lambda m: f'url(#{prefix}{m.group(1)})', svg_text)
    svg_text = re.sub(r'(xlink:href|href)="#([^"]+)"', lambda m: f'{m.group(1)}="#{prefix}{m.group(2)}"', svg_text)
    return svg_text


def _svg_size(root_tag):
    def attr(name):
        m = re.search(rf'\b{name}="([^"]+)"', root_tag)
        return m.group(1) if m else None
    vb = attr('viewBox')
    if vb:
        vx, vy, vw, vh = (float(v) for v in vb.replace(',', ' ').split())
    else:
        vx = vy = 0.0
        vw = float(re.sub(r'[a-z]+$', '', attr('width')))
        vh = float(re.sub(r'[a-z]+$', '', attr('height')))
    return vx, vy, vw, vh


def nest_svg(path, x, y, w, h, prefix, align='xMinYMin'):
    """Return the panel SVG at `path` as a nested <svg> element filling box (x, y, w, h)."""
    text = _prefix_ids(_strip_prolog(Path(path).read_text()), prefix)
    m = re.search(r'<svg\b[^>]*>', text)
    root = m.group(0)
    vx, vy, vw, vh = _svg_size(root)
    body = text[m.end():text.rfind('</svg>')]
    # Keep the panel's namespace declarations; drop its size/position attributes
    ns = ' '.join(re.findall(r'xmlns(?::\w+)?="[^"]+"', root))
    return (f'<svg {ns} x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" '
            f'viewBox="{vx} {vy} {vw} {vh}" preserveAspectRatio="{align} meet" overflow="visible">'
            f'{body}</svg>')


def text_el(s, x, y, size, font, weight='normal', anchor='start', fill='#000000', rotate=None):
    tr = f' transform="rotate({rotate} {x:.3f} {y:.3f})"' if rotate else ''
    return (f'<text x="{x:.3f}" y="{y:.3f}" font-family="{font}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}"{tr}>{escape(s)}</text>')


def line_el(x1, y1, x2, y2, stroke='#000000', width=0.5):
    return (f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def compose(items, width_pt, height_pt, out_svg, scale=1.0):
    """Write the page SVG. Each item is an SVG fragment string.

    The items are in page points; the page is written `scale` times that size (the viewBox
    keeps the page-point coordinates).
    """
    head = (f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{round(width_pt * scale, 3)}pt" height="{round(height_pt * scale, 3)}pt" '
            f'viewBox="0 0 {width_pt} {height_pt}">\n'
            f'<rect x="0" y="0" width="{width_pt}" height="{height_pt}" fill="#FFFFFF"/>\n')
    Path(out_svg).write_text(head + '\n'.join(items) + '\n</svg>\n')


def render(svg_path, pdf_path=None, png_path=None, png_size=None):
    """Render the page SVG to PDF and PNG with cairosvg (fonts resolved through fontconfig)."""
    import cairosvg
    if pdf_path:
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
    if png_path:
        kw = {}
        if png_size:
            kw = {'output_width': png_size[0], 'output_height': png_size[1]}
        cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), background_color='white', **kw)
