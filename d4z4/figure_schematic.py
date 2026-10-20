#!/usr/bin/env python3
"""
Fig. S9A: the classification-pipeline schematic, as published.

This is the layout of create_pipeline_diagram.py as it was reworked for the manuscript: the
seven steps on the left, three worked examples on the right, the four classification rules
as a 2 x 2 grid and a small legend. Text and geometry follow the published figure;
colours are those of create_pipeline_diagram.py (which match plot_d4z4_arrays.py).

Coordinates are in points on the published page (as read from the Illustrator file) and are
shifted by `origin` when drawn.
"""

from xml.sax.saxutils import escape

COLORS = {
    'process': '#FFF3E0', 'edge': '#666666', 'arrow': '#455A64', 'text': '#000000',
    'bsat': '#fa99ff', 'exon': '#090088', 'intron': '#198450', 'intergenic': '#E89EB8',
    'ct': '#B8AFA2', 'plam': '#228B22', 'line': '#FFB366', 'sine': '#66B3FF',
    'canonical': '#90EE90', 'canonical_ct': '#87CEEB', 'degraded': '#FF6B6B', 'transition': '#FFD700',
}

SIZE = 9.31          # all schematic text
LINE_PITCH = 11.31   # baseline-to-baseline within a box

STEPS = [
    ('Step 1: Extract β-satellite array', ['bsat_specific calls from KaryoScope']),
    ('Step 2: Merge Clusters', ['Merge fragments (gap ≤ 400bp)']),
    ('Step 3: Annotate Clusters', ['Calculate exon fraction from gene track']),
    ('Step 4: Identify Arrays', ['Group exon-dense clusters', '(≥50% exon, ≥1500bp)']),
    ('Step 5: Extend Arrays', ['Add trailing fragments (gap ≤ 1500bp)']),
    ('Step 6: Find Terminals', ['Search for terminal β-satellite', 'or ct_specific terminal']),
    ('Step 7: Classify', ['Apply classification rules']),
]
STEP_X = (30.86, 195.97)
STEP_TOP0, STEP_PITCH, STEP_H = 5.67, 62.23, 35.68

RULES = [  # title, lines, colour, (x0, y0)
    ('Canonical (C+)', ['β-sat terminal + ≥1 elem'], 'canonical', (227.4, 310.5)),
    ('Canonical (C-)', ['ct_specific terminal', '+ ≥3 elem'], 'canonical_ct', (356.15, 310.5)),
    ('Degraded', ['≥3 elements + no terminal'], 'degraded', (227.4, 374.31)),
    ('Transition element', ['isolated elements'], 'transition', (356.15, 374.31)),
]
RULE_W, RULE_H = 115.5, 39.0

# Worked examples: D4Z4 unit k starts at UNIT_X0 + k * UNIT_PITCH
UNIT_X0, UNIT_PITCH = 262.01, 40.67
ROW_LABEL_RIGHT = 253.4
ROW_H = 10.37
CAP_HEIGHT = 0.716   # Arial/Helvetica cap (and digit) height, in em


def _centred_baseline(top, height):
    """Baseline that centres capitals and digits vertically in a box."""
    return top + height / 2 + CAP_HEIGHT * SIZE / 2


def _rect(x, y, w, h, fill, stroke=None, sw=0, rx=0):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ''
    r = f' rx="{rx}"' if rx else ''
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}"{s}{r}/>'


def _text(s, x, y, font, bold=False, anchor='start', fill=None):
    w = ' font-weight="bold"' if bold else ''
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-family="{font}" font-size="{SIZE}"{w} '
            f'text-anchor="{anchor}" fill="{fill or COLORS["text"]}">{escape(s)}</text>')


def _unit_satellite(k, n=3):
    """The bsat_specific fragments of one D4Z4 unit (three 8.13 pt blocks at 12.2 pt pitch)."""
    x = UNIT_X0 + k * UNIT_PITCH
    return [(x + i * 12.2, 8.13, 'bsat') for i in range(n)]


def _unit_gene(k, partial=False):
    """Gene track of one unit: exon, intron, exon (the exon-dense D4Z4 ORF)."""
    x = UNIT_X0 + k * UNIT_PITCH
    segs = [(x, 16.27, 'exon'), (x + 16.27, 8.13, 'intron'), (x + 24.4, 12.2, 'exon')]
    return segs[:2] if partial else segs


def _example(out, font, title, title_y, rows, array, extra):
    """rows: {label: (y_top, [(x, w, colour), ...])}; array: (y_top, n_units, box colour, box x1)."""
    out.append(_text(title, 257.95, title_y, font, bold=True))
    for label, (y, segs) in rows.items():
        out.append(_text(label, ROW_LABEL_RIGHT, y + 6.48, font, anchor='end'))
        for x, w, c in segs:
            out.append(_rect(x, y, w, ROW_H, COLORS[c]))
    y, n, color, x1 = array
    out.append(_text('D4Z4 arrays', ROW_LABEL_RIGHT, y + 10.42, font, anchor='end'))
    out.append(_rect(259.17, y, x1 - 259.17, 18.25, COLORS[color], COLORS['edge'], 2.0, rx=3))
    for k in range(n):
        x = UNIT_X0 + k * UNIT_PITCH
        out.append(_rect(x, y + 4.15, 32.54, 9.95, COLORS['bsat'], COLORS['edge'], 0.5))
        out.append(_text(str(k + 1), x + 16.27, _centred_baseline(y + 4.15, 9.95), font, anchor='middle'))
    out.extend(extra)


def svg_fragment(origin, font, size=None):
    """The schematic as an SVG <g>, shifted so that page point `origin` maps to (0, 0).

    size: text size in page points (default SIZE, the published 9.31).
    """
    global SIZE
    SIZE = size or 9.31
    ox, oy = origin
    out = [f'<g transform="translate({-ox:.2f} {-oy:.2f})">']

    # Left: the seven steps, joined by arrows
    x0, x1 = STEP_X
    cx = (x0 + x1) / 2
    for i, (title, lines) in enumerate(STEPS):
        top = STEP_TOP0 + i * STEP_PITCH
        out.append(_rect(x0, top, x1 - x0, STEP_H, COLORS['process'], COLORS['edge'], 1.5, rx=4))
        first = top + 14.62 - (len(lines) - 1) * 4.12
        out.append(_text(title, cx, first, font, bold=True, anchor='middle'))
        for j, line in enumerate(lines, 1):
            out.append(_text(line, cx, first + j * LINE_PITCH, font, anchor='middle'))
        if i < len(STEPS) - 1:
            y_a, y_b = top + STEP_H + 2.0, top + STEP_PITCH - 3.0
            out.append(f'<path d="M{cx:.2f},{y_a:.2f} L{cx:.2f},{y_b:.2f} M{cx - 1.15:.2f},{y_b - 2.36:.2f} '
                       f'L{cx:.2f},{y_b:.2f} L{cx + 1.15:.2f},{y_b - 2.36:.2f}" fill="none" '
                       f'stroke="{COLORS["arrow"]}" stroke-width="1.5"/>')

    # Example 1: canonical C+, beta-satellite terminal
    sat = [s for k in range(5) for s in _unit_satellite(k)] + _unit_satellite(5, 2) + \
          [(501.95 + i * 12.2, 8.13, 'bsat') for i in range(3)]
    gene = [s for k in range(5) for s in _unit_gene(k)] + \
           [(493.82, 8.13, 'intron'), (501.95, 8.13, 'intergenic'), (510.08, 24.4, 'intron')]
    rep = [(469.42, 12.2, 'line'), (542.62, 12.2, 'sine'), (562.95, 12.2, 'line')]
    _example(out, font, 'Example 1: Canonical (C+) - β-satellite terminal', 16.54,
             {'Satellites': (26.63, sat), 'Genes': (43.22, gene), 'Repeats': (59.81, rep)},
             (72.47, 5, 'canonical', 468.2),
             [_text('intron-rich', 515.59, 61.88, font, anchor='middle'),
              _rect(493.82, 76.62, 40.66, 9.95, COLORS['plam'], COLORS['edge'], 0.5),
              _text('β-sat', 514.15, _centred_baseline(76.62, 9.95), font, anchor='middle')])

    # Example 2: canonical C-, ct_specific terminal
    sat = [s for k in range(4) for s in _unit_satellite(k)] + _unit_satellite(4, 2) + \
          [(481.62, 101.67, 'ct')]
    gene = [s for k in range(4) for s in _unit_gene(k)] + _unit_gene(4, partial=True) + \
           [(481.62, 20.33, 'exon'), (501.95, 20.33, 'intergenic'), (522.28, 20.34, 'exon'),
            (542.62, 20.33, 'intergenic')]
    rep = [(453.15, 16.27, 'line'), (583.29, 20.33, 'sine')]
    _example(out, font, 'Example 2: Canonical (C-) - ct_specific terminal', 119.65,
             {'Satellites': (130.34, sat), 'Genes': (146.93, gene), 'Repeats': (163.53, rep)},
             (176.18, 4, 'canonical_ct', 447.86),
             [_text('ct_specific', 533.92, 138.22, font, anchor='middle')])

    # Example 3: degraded, no terminal
    sat = [s for k in range(3) for s in _unit_satellite(k)]
    gene = [s for k in range(3) for s in _unit_gene(k)] + [(400.28, 61.0, 'intergenic')]
    _example(out, font, 'Example 3: Degraded (D) - no terminal detected', 224.68,
             {'Satellites': (233.55, sat), 'Genes': (250.14, gene)},
             (262.79, 3, 'degraded', 382.8),
             [_text('No terminal', 430.48, 246.23, font, anchor='middle')])

    # Classification rules, 2 x 2
    for title, lines, color, (x, y) in RULES:
        out.append(_rect(x, y, RULE_W, RULE_H, COLORS[color], COLORS['edge'], 1.0, rx=4))
        first = y + 16.75 - (len(lines) - 1) * 4.13
        out.append(_text(title, x + RULE_W / 2, first, font, bold=True, anchor='middle'))
        for j, line in enumerate(lines, 1):
            out.append(_text(line, x + RULE_W / 2, first + j * LINE_PITCH, font, anchor='middle'))

    # Legend
    legend = [('Satellites', 216.55, [('bsat_specific', 'bsat', 231.33)]),
              ('Genes', 253.17, [('exon', 'exon', 265.17), ('intron', 'intron', 282.06),
                                 ('intergenic', 'intergenic', 299.84)]),
              ('Terminal region', 318.9, [('ct_specific terminal', 'ct', 331.92),
                                          ('β-sat terminal', 'plam', 350.26)]),
              ('Repeats', 375.85, [('LINE', 'line', 387.64), ('SINE', 'sine', 406.0)])]
    for header, hy, items in legend:
        out.append(_text(header, 493.82, hy, font, bold=True))
        for label, c, y in items:
            out.append(_rect(494.5, y - 7.39, 16.27, 9.95, COLORS[c], COLORS['edge'], 0.5))
            out.append(_text(label, 516.7, y, font))
    out.append('</g>')
    return '\n'.join(out)
