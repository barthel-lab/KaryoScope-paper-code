#!/usr/bin/env python3
"""
chr4 D4Z4 sunburst (HPRC2 Supplementary Fig. D4Z4_7f, formerly Fig. 7f), with and without
the assembly-QC filter (see make_figures.py --figures sunburst).

Counts come from plot_array_distribution.calculate_distributions(), the function behind the
published panel, so the rings are those of plot_array_distribution.py:
    inner   number of arrays (0, 1, 2+)
    second  single-array haplotypes: 4qA (beta-satellite terminal) or 4qB
    third   4qA: pLAM functional (ATTAAA), non-functional (ATCAAA) or absent
    outer   functional pLAM: repeat units 1-10, 11-20, 21-30, 31+
The ring radii and colours are those of plot_array_distribution.plot_sunburst().

The labels are laid out by rule. layout() computes every wedge and label from the counts and a
few parameters (start angle, direction, order of each ring, how much a label written along its
ring may tilt, whether labels at the top and bottom sit centred above/below the chart); score()
measures a layout (mean tilt of the label text from horizontal, area of the rectangle around
chart, labels and title) and rejects any in which a label touches a ring or another label;
render() draws it. A label is written along its ring when it fits in its wedge (shrinking to
MIN_LABEL_PT at most; it slides along its wedge to the most upright position) and tilts no
more than max_tilt (or is listed in inside, which keeps it in
its ring regardless); otherwise it is written horizontally outside, with a leader line, and
outside labels on one side are spread so they never overlap. label_style chooses the text: name
and count (full), count only (values), or count in the rings with the name on leader labels,
whose wedges are too thin to identify by colour (hybrid). Neighbouring outside labels from
different rings are set GROUP_GAP further apart, so each ring's labels read as a group.
The parameters in figures_config.yaml (sunburst:) were chosen by searching every combination
with score() (search_layouts()).

Writes sunburst_chr4.{svg,pdf,png} (all haplotypes), sunburst_chr4_qc.{svg,pdf,png} (only
haplotypes passing the array-focused assembly QC, qc_pass in the classification table) and
sunburst_chr4_counts.tsv (every count drawn, both versions).
"""

import functools
import itertools
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle, Wedge
from matplotlib.textpath import TextToPath

import plot_array_distribution as pad

# Geometry and colours of plot_array_distribution.plot_sunburst()
RING_WIDTH = 0.19
RADII = {1: 0.25, 2: 0.47, 3: 0.69, 4: 0.91}   # inner radius of each ring
COLORS = {
    'No arrays': '#CCCCCC', '1 array': '#98FB98', '2+ arrays': '#DDA0DD',
    'C+': '#90EE90', 'C-': '#87CEEB', 'D': '#FF6B6B',
    'Functional': '#40D392', 'Non-functional': '#F07167', 'Absent': '#6B7280',
    '1-10 (FSHD1)': '#FF6B6B', '11-20 (FSHD2)': '#FFA07A', '21-30': '#98FB98', '31+': '#87CEEB',
}
# palette 'ordinal': the two ordered rings as one-hue ramps, darker for more (array count:
# amber, 0 lightest; repeat units: violet, 31+ darkest), in hues no other ring uses, so no
# colour repeats across rings. The 4qA/4qB and pLAM colours stay those of the D4Z4 track views
# (A-like, B-like, pLAM marks).
ORDINAL_COLORS = dict(COLORS, **{
    'No arrays': '#FBF1D3', '1 array': '#F4D78A', '2+ arrays': '#C98A00',
    '1-10 (FSHD1)': '#DCD7F7', '11-20 (FSHD2)': '#B2A8EE', '21-30': '#8676E2', '31+': '#5441C2',
})
PALETTES = {'published': COLORS, 'ordinal': ORDINAL_COLORS}
# Display names (the counting function's keys stay as they are)
NAMES = {'C+': '4qA', 'C-': '4qB', 'D': 'Degraded', '1-10 (FSHD1)': '1–10', '11-20 (FSHD2)': '11–20',
         '21-30': '21–30', '31+': '31+'}
LEGEND = [
    ('Inner ring, array count:', [('No arrays', 'No arrays'), ('1 array', '1 array'), ('2+ arrays', '2+ arrays')]),
    ('Second ring, array type (1 array only):', [('C+', '4qA (with β-satellite)'), ('C-', '4qB (no β-satellite)'),
                                                 ('D', 'Degraded')]),
    ('Third ring, pLAM status (4qA only):', [('Functional', 'Functional pLAM'), ('Non-functional', 'Non-functional pLAM'),
                                            ('Absent', 'Absent')]),
    ('Outer ring, repeat count (functional pLAM only):', [('1-10 (FSHD1)', '1–10 (FSHD1 range)'),
                                                          ('11-20 (FSHD2)', '11–20 (intermediate)'),
                                                          ('21-30', '21–30 (normal)'), ('31+', '31+ (long)')]),
]
LEVELS = {1: ('level1', ('No arrays', '1 array', '2+ arrays')),
          2: ('level2', ('C+', 'C-', 'D')),
          3: ('level3', ('Functional', 'Non-functional', 'Absent')),
          4: ('level4', ('1-10 (FSHD1)', '11-20 (FSHD2)', '21-30', '31+'))}
PARENT = {2: '1 array', 3: 'C+', 4: 'Functional'}   # the wedge each ring subdivides

DEFAULT_PARAMS = {'start': 90.0, 'clockwise': True, 'max_tilt': 90.0, 'polar_labels': False,
                  'label_style': 'full', 'inside': (), 'label_arrangement': 'column', 'center_gap': False,
                  'other_step': 1.0, 'outer_spacing': 'angle',
                  'order1': LEVELS[1][1], 'order2': LEVELS[2][1], 'order3': LEVELS[3][1], 'order4': LEVELS[4][1]}

PT_PER_UNIT = 80.0              # points per data unit (rings 0.25-1.10 units: 2.4 in across)
UNIT = 1 / PT_PER_UNIT          # data units per point
LABEL_PT, LEGEND_PT, TITLE_PT = 7.5, 8.0, 9.0
MIN_LABEL_PT = 6.5             # a ring label may shrink this far to stay inside its wedge
EDGE_PT = 1.0                   # wedge outline
CAP, DESC = 0.716, 0.21         # Arial cap height and descent / font size
LINE_H = 1.2 * LABEL_PT * UNIT  # line height of an outside label
LABEL_GAP = 0.09                # outside labels: distance beyond the local rings
LEADER_MARGIN = 0.04            # a leader's second line keeps this far from every ring
GROUP_GAP = 0.08                # extra space between neighbouring outside labels from different rings
FONT = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']
_t2p = TextToPath()


@functools.lru_cache(maxsize=None)
def text_width_pt(s, size, bold=False):
    fp = FontProperties(family=FONT, size=size, weight='bold' if bold else 'normal')
    return _t2p.get_text_width_height_descent(s, fp, ismath=False)[0]


def pct_text(value, total):
    p = 100 * value / total if total else 0
    return '<1%' if 0 < p < 0.5 else f'{p:.0f}%'


def in_span(theta, a, b):
    """theta within the wedge from a to b (degrees, either direction)."""
    lo, hi = min(a, b), max(a, b)
    t = (theta - lo) % 360
    return t <= hi - lo


# -- layout -------------------------------------------------------------------------------------

def segments(counts, order, start, span, sign):
    """[(key, value, a, b)] from start over span degrees (sign -1 clockwise), zero counts skipped."""
    total = sum(counts.get(k, 0) for k in order)
    out, a = [], start
    for k in order:
        v = counts.get(k, 0)
        if v == 0 or total == 0:
            continue
        b = a + sign * span * v / total
        out.append((k, v, a, b))
        a = b
    return out


def arc_chars(s, ring, mid, size):
    """Per-character (ch, width, x, y, rotation) for s written along a ring, centred on mid, in
    whichever reading direction (letter tops outward or inward) keeps the text more upright."""
    a = _arc_chars(s, ring, mid, size, True)
    b = _arc_chars(s, ring, mid, size, False)
    return min(a, b, key=mean_tilt)


def mean_tilt(chars):
    return sum(w * tilt(r) for _, w, _, _, r in chars) / sum(w for _, w, _, _, _ in chars)


def _arc_chars(s, ring, mid, size, clockwise):
    r = RADII[ring] + RING_WIDTH / 2
    widths = [text_width_pt(ch, size) * UNIT for ch in s]
    span = sum(widths) / r
    cap = CAP * size * UNIT
    rb = r - cap / 2 if clockwise else r + cap / 2
    theta = math.radians(mid) + (span / 2 if clockwise else -span / 2)
    out = []
    for ch, w in zip(s, widths):
        step = w / rb
        c = theta - step / 2 if clockwise else theta + step / 2
        rot = math.degrees(c) - 90 if clockwise else math.degrees(c) + 90
        out.append((ch, w, rb * math.cos(c), rb * math.sin(c), rot))
        theta = theta - step if clockwise else theta + step
    return out


def best_position(s, ring, a, b, size, step=2.0):
    """s along its ring, centred on the vertical axis (90 or 270 degrees) if it fits there, else at
    the most upright position that keeps it inside the wedge a-b (ties go to the position nearest
    the wedge's middle): a label need not sit at the middle of its wedge."""
    r = RADII[ring] + RING_WIDTH / 2
    half = math.degrees((text_width_pt(s, size) * UNIT + 0.16) / 2 / r)   # stop short of the wedge edges
    lo, hi = min(a, b) + half, max(a, b) - half
    mid = (a + b) / 2
    # centred on the vertical line through the centre (top or bottom) whenever the label fits there;
    # nearest the wedge's middle, the bottom on a tie
    axis = sorted((abs(c - mid), 0 if c % 360 == 270 else 1, c) for base in (90, 270)
                  for c in (base + 360 * k for k in range(-2, 3)) if lo - 1e-9 <= c <= hi + 1e-9)
    if axis:
        return arc_chars(s, ring, axis[0][2], size)
    best = None
    c = lo
    while c <= hi + 1e-9:
        chars = arc_chars(s, ring, c, size)
        key = (round(mean_tilt(chars), 1), abs(c - mid))
        if best is None or key < best[0]:
            best = (key, chars)
        c += step
    return best[1] if best else arc_chars(s, ring, mid, size)


def tilt(rot):
    """Angle of a letter from upright, 0-180 degrees (upside-down text counts as 180)."""
    return abs((rot + 180) % 360 - 180)


def fit_size(s, ring, a, b):
    r = RADII[ring] + RING_WIDTH / 2
    arc = math.radians(abs(a - b)) * r
    size = LABEL_PT
    while size >= MIN_LABEL_PT - 1e-9:
        if text_width_pt(s, size) * UNIT + 0.08 <= arc:
            return size
        size -= 0.25
    return None


def layout(data, title, params=None):
    """Every wedge and label of the chart (data units), or None if no wedge exists."""
    p = dict(DEFAULT_PARAMS, **(params or {}))
    sign = -1 if p['clockwise'] else 1
    if p.get('center_gap'):
        # rotate so the part of the second ring left empty (the inner ring outside '1 array') is
        # centred on the positive x axis
        one = [(a, b) for k, _, a, b in segments(data['level1'], p['order1'], p['start'], 360.0, sign)
               if k == '1 array']
        if one:
            gap_mid = (one[0][0] + one[0][1]) / 2 + 180
            p['start'] -= (gap_mid + 180) % 360 - 180
    L = {'wedges': [], 'arc': [], 'outside': [], 'title': title, 'center': f'n = {data["total"]}'}
    spans = {}
    parent_total = {1: data['total']}
    for ring in (1, 2, 3, 4):
        key, _ = LEVELS[ring]
        counts = data[key]
        if ring == 1:
            start, span = p['start'], 360.0
        else:
            if PARENT[ring] not in spans:
                break
            a0, b0 = spans[PARENT[ring]]
            start, span = a0, abs(b0 - a0)
        n = sum(counts.values())
        parent_total[ring] = n
        for k, v, a, b in segments(counts, p[f'order{ring}'], start, span, sign):
            spans[k] = (a, b)
            L['wedges'].append((ring, a, b, k))
            name = NAMES.get(k, k)
            pct = pct_text(v, n)
            if p['label_style'] == 'values':          # counts only; the legend names the wedges
                inside = outside = f'n = {v} ({pct})'
            elif p['label_style'] == 'hybrid':        # counts in the rings, names on the leader labels
                inside = f'n = {v} ({pct})'
                outside = f'{name}: n = {v} ({pct})' if ring in (1, 4) else f'{name}\nn = {v} ({pct})'
            else:
                inside = f'{name}: n = {v} ({pct})'
                outside = inside if ring in (1, 4) else f'{name}\nn = {v} ({pct})'
            size = fit_size(inside, ring, a, b)
            must = k in p['inside']                   # this wedge's label must stay in its ring
            if size:
                chars = best_position(inside, ring, a, b, size)
                t = mean_tilt(chars)
                if must or t <= p['max_tilt'] + 1e-9:
                    L['arc'].append({'text': inside, 'size': size, 'chars': chars, 'tilt': t})
                    continue
            if must:
                L['invalid'] = True
            L['outside'].append({'mid': (a + b) / 2, 'r_edge': RADII[ring] + RING_WIDTH, 'text': outside,
                                 'ring': ring})
    L['other_step'] = p['other_step']
    L['outer_spacing'] = p['outer_spacing']
    place_outside(L, p['polar_labels'], p['label_arrangement'])
    return L


def outer_extent(L, theta, window=8.0):
    """Outer radius of the rings drawn within +-window degrees of theta."""
    rs = [RADII[1] + RING_WIDTH]
    for ring, a, b, _ in L['wedges']:
        if any(in_span(theta + d, a, b) for d in (-window, 0, window)):
            rs.append(RADII[ring] + RING_WIDTH)
    return max(rs)


def text_box(lab):
    """(x0, y0, x1, y1) of an outside label."""
    lines = lab['text'].split('\n')
    w = max(text_width_pt(s, LABEL_PT) for s in lines) * UNIT
    h = len(lines) * LINE_H
    x, y = lab['x'], lab['y']
    va = lab.get('va', 'center')
    if va == 'bottom':                    # label sits on its leader end
        y = y + h / 2
    elif va == 'top':                     # label hangs from its leader end
        y = y - h / 2
    if lab['ha'] == 'left':
        x0 = x
    elif lab['ha'] == 'right':
        x0 = x - w
    else:
        x0 = x - w / 2
    return (x0, y - h / 2, x0 + w, y + h / 2)


def place_on_circle(items, side, f=1.0):
    """Side labels whose leader lines end on one circle around the centre (the text starts just
    beyond the end point; each leader runs out from its wedge, then across). The outermost ring's
    labels split that ring's arc into equal parts, one end point at the middle of each. Each other
    ring's labels (a
    group) are evenly spaced in angle, the step fitted to the angles of the wedges they label
    (least squares) scaled by f between the tightest spacing (0) and the full fit (1), widened
    only as needed to keep neighbouring labels clear; labels from different rings are GROUP_GAP
    further apart. False if no arrangement fits."""
    sgn = 1 if side == 'right' else -1
    R = max(lab['r'] for lab in items)
    gap = math.degrees(GROUP_GAP / R)
    # angles from the horizontal on this side, positive up
    mids = [math.degrees(math.atan2(math.sin(math.radians(lab['mid'])),
                                    sgn * math.cos(math.radians(lab['mid'])))) for lab in items]
    groups, start = [], 0
    for i in range(1, len(items) + 1):
        if i == len(items) or items[i]['ring'] != items[i - 1]['ring']:
            groups.append(list(range(start, i)))
            start = i
    steps, t0s = [], []
    for g in groups:
        n = len(g)
        if n > 1:                                   # fit mids ~ t0 - step * k
            ks = list(range(n))
            km, mm = sum(ks) / n, sum(mids[i] for i in g) / n
            slope = -sum((k - km) * (mids[i] - mm) for k, i in zip(ks, g)) / sum((k - km) ** 2 for k in ks)
            step = 2.0 + f * max(0.0, slope - 2.0)
        else:
            step = 2.0
        steps.append(step)
        t0s.append(sum(mids[i] + step * k for k, i in enumerate(g)) / n)
    # the outermost ring's labels split that ring's arc evenly: each end point at the middle of
    # one of n equal parts, top to bottom
    outer = max(r for r, _, _, _ in L_CURRENT[0]['wedges'])
    edges = [math.degrees(math.atan2(math.sin(math.radians(t)), sgn * math.cos(math.radians(t))))
             for r, a, b, _ in L_CURRENT[0]['wedges'] if r == outer for t in (a, b)]
    for gi, g in enumerate(groups):
        if items[g[0]]['ring'] == outer and edges:
            top, bottom = max(edges), min(edges)
            steps[gi] = (top - bottom) / len(g)
            t0s[gi] = top - steps[gi] / 2
    # the other rings' labels continue along the circle at other_step times the outer ring's step,
    # one GROUP_GAP further apart where the ring changes
    oi = next((gi for gi, g in enumerate(groups) if items[g[0]]['ring'] == outer), None)
    if oi is not None:
        st = steps[oi] * L_CURRENT[0].get('other_step', 1.0)
        for gi in range(oi + 1, len(groups)):            # below the outer ring's labels
            steps[gi] = st
            last = t0s[gi - 1] - steps[gi - 1] * (len(groups[gi - 1]) - 1)
            t0s[gi] = last - st - gap
        for gi in range(oi - 1, -1, -1):                 # above them
            steps[gi] = st
            t0s[gi] = t0s[gi + 1] + st + gap + st * (len(groups[gi]) - 1)

    vertical = L_CURRENT[0].get('outer_spacing') == 'vertical'
    centre = {}                                     # vertical: each outer label's target centre y

    def angles():
        ts = [0.0] * len(items)
        for g, st, t0 in zip(groups, steps, t0s):
            for k, i in enumerate(g):
                ts[i] = t0 - st * k
            if vertical and items[g[0]]['ring'] == outer and len(g) > 1:
                # same top and bottom positions, the labels' centres equally spaced in height
                y0, y1 = R * math.sin(math.radians(ts[g[0]])), R * math.sin(math.radians(ts[g[-1]]))
                for k, i in enumerate(g):
                    yc = y0 + (y1 - y0) * k / (len(g) - 1)
                    centre[id(items[i])] = yc
                    ts[i] = math.degrees(math.asin(max(-1.0, min(1.0, yc / R))))
        return ts

    def need(i):                                    # vertical clearance between items i and i+1
        extra = GROUP_GAP if items[i]['ring'] != items[i + 1]['ring'] else 0.0
        return (items[i]['h'] + items[i + 1]['h']) / 2 + 0.01 + extra   # same padding as _overlap()

    group_of = {i: gi for gi, g in enumerate(groups) for i in g}
    for _ in range(2000):
        ts = angles()
        if max(ts) > 88 or min(ts) < -88:
            return False
        bad = None
        for i in range(len(items) - 1):
            dy = R * (math.sin(math.radians(ts[i])) - math.sin(math.radians(ts[i + 1])))
            same = group_of[i] == group_of[i + 1]
            if dy < need(i):
                bad = (i, same)
                break
        if bad is None:
            break
        i, same = bad
        if same:
            steps[group_of[i]] += 0.25              # widen this group's step about its centre
            g = groups[group_of[i]]
            t0s[group_of[i]] += 0.25 * (len(g) - 1) / 2
        else:                                       # move the lower groups down
            for gi in range(group_of[i + 1], len(groups)):
                t0s[gi] -= 0.25
    else:
        return False
    bend = {}                                       # bend radius each label needs on its own
    for lab, t in zip(items, angles()):
        a, tr = math.radians(lab['mid']), math.radians(t)
        xe, ye = sgn * R * math.cos(tr), R * math.sin(tr)               # leader end, on the circle
        r_out = outer_extent(L_CURRENT[0], lab['mid']) + 0.03
        lab['x'], lab['y'], lab['ha'] = xe + 0.02 * sgn, ye, 'left' if side == 'right' else 'right'
        lab.pop('_va', None)
        if id(lab) in centre:
            # keep the label's centre at its target height: centred on its end point if the leader
            # clears the label, else sitting on (or hanging from) an end point half a line lower
            # (higher), still on the circle
            yc = centre[id(lab)]
            for va, off in (('center', 0.0), ('bottom', lab['h'] / 2), ('top', -lab['h'] / 2)):
                ye = yc - off
                xe = sgn * math.sqrt(max(R * R - ye * ye, 0.0))
                lab['x'], lab['y'], lab['va'] = xe + 0.02 * sgn, ye, va
                lab['leader'] = [((lab['r_edge'] + 0.015) * math.cos(a), (lab['r_edge'] + 0.015) * math.sin(a)),
                                 (r_out * math.cos(a), r_out * math.sin(a)), (xe, ye)]
                if not _crosses_own_label(lab):
                    break
            lab['_va'] = va
        lab['_end'] = (xe, ye)
        # two lines: out from the wedge to the bend, then across to the end point; the bend is the
        # nearest point on the wedge's radial line from which both lines clear every ring
        p0 = ((lab['r_edge'] + 0.015) * math.cos(a), (lab['r_edge'] + 0.015) * math.sin(a))
        rb = lab['r_edge'] + 0.04
        while rb < r_out:
            e = (rb * math.cos(a), rb * math.sin(a))
            if not any(_in_ring(p0[0] + (e[0] - p0[0]) * k / 40, p0[1] + (e[1] - p0[1]) * k / 40, L_CURRENT[0])
                       for k in range(1, 41)) and \
                    not any(_in_ring(e[0] + (xe - e[0]) * k / 40, e[1] + (ye - e[1]) * k / 40, L_CURRENT[0],
                                     margin=LEADER_MARGIN) for k in range(1, 41)):
                break
            rb += 0.01
        bend[id(lab)] = min(rb, r_out)
    for g in groups:                                # one bend radius per ring, so its leaders run parallel
        rg = max(bend[id(items[i])] for i in g)
        for i in g:
            lab = items[i]
            a = math.radians(lab['mid'])
            xe, ye = lab.pop('_end')
            elbow = (rg * math.cos(a), rg * math.sin(a))
            lab['leader'] = [((lab['r_edge'] + 0.015) * math.cos(a), (lab['r_edge'] + 0.015) * math.sin(a)),
                             elbow, (xe, ye)]
            # the label is centred on its end point unless its leader would cut through it; then it
            # sits on the end point (leader rising into it) or hangs from it (leader falling into it)
            if '_va' in lab:                          # chosen above to keep its centre
                lab['va'] = lab.pop('_va')
                continue
            lab['va'] = 'center'
            if _crosses_own_label(lab):
                lab['va'] = 'bottom' if ye > elbow[1] else 'top'
    return True


L_CURRENT = [None]   # the layout being placed (outer_extent needs its wedges)


def _crosses_own_label(lab):
    """True if the leader passes through its own label anywhere but at its end point."""
    b = text_box(lab)
    segs = list(zip(lab['leader'], lab['leader'][1:]))
    end = segs[-1][1]
    for p, q in segs:
        for k in range(21):
            x, y = p[0] + (q[0] - p[0]) * k / 20, p[1] + (q[1] - p[1]) * k / 20
            if math.hypot(x - end[0], y - end[1]) > 0.04 and b[0] <= x <= b[2] and b[1] <= y <= b[3]:
                return True
    return False


def place_outside(L, polar, arrangement='column'):
    """Positions of the outside labels and their leader lines. arrangement 'column': each side's
    labels stacked beside their wedges, spread just enough not to overlap; 'circle': their start
    points on one circle at even angular steps (place_on_circle())."""
    L_CURRENT[0] = L
    labels = L['outside']
    for lab in labels:
        a = math.radians(lab['mid'])
        lab['r'] = outer_extent(L, lab['mid']) + LABEL_GAP
        lab['h'] = (lab['text'].count('\n') + 1) * LINE_H
        c = math.cos(a)
        lab['mode'] = 'polar' if polar and abs(c) < 0.26 else ('right' if c >= 0 else 'left')
    for side in ('right', 'left'):
        items = sorted([lab for lab in labels if lab['mode'] == side],
                       key=lambda lab: -lab['r'] * math.sin(math.radians(lab['mid'])))
        if arrangement == 'circle' and items:
            best = None                             # the even spacing with the smallest chart
            for f10 in range(11):
                if place_on_circle(items, side, f10 / 10):
                    sc = score(L, partial=True)
                    if sc and (best is None or sc[0] < best[0] - 1e-9):
                        best = (sc[0], f10 / 10)
            if best is not None:
                place_on_circle(items, side, best[1])
                continue
        for lab in items:
            lab['y'] = lab['r'] * math.sin(math.radians(lab['mid']))
        for _ in range(200):                                     # relax overlaps, moving both apart
            moved = False
            for up, dn in zip(items, items[1:]):
                need = (up['h'] + dn['h']) / 2 + 0.02 + (GROUP_GAP if up['ring'] != dn['ring'] else 0)
                gap = up['y'] - dn['y']
                if gap < need:
                    d = (need - gap) / 2
                    up['y'] += d
                    dn['y'] -= d
                    moved = True
            if not moved:
                break
        sgn = 1 if side == 'right' else -1
        for lab in items:
            a, y, r = math.radians(lab['mid']), lab['y'], lab['r']
            x = sgn * math.sqrt(max(r * r - y * y, 0.0))         # on (or above) the label circle
            lab['x'], lab['ha'] = x + 0.04 * sgn, 'left' if side == 'right' else 'right'   # text extends outward
            lab['leader'] = [((lab['r_edge'] + 0.015) * math.cos(a), (lab['r_edge'] + 0.015) * math.sin(a)),
                             ((r - 0.05) * math.cos(a), (r - 0.05) * math.sin(a)),
                             (x + 0.02 * sgn, y)]
    for lab in labels:
        if lab['mode'] != 'polar':
            continue
        a, r = math.radians(lab['mid']), lab['r']
        up = math.sin(a) >= 0
        x = r * math.cos(a)
        y = r * math.sin(a) + (lab['h'] / 2 + 0.02) * (1 if up else -1)
        lab['x'], lab['y'], lab['ha'] = x, y, 'center'
        lab['leader'] = [((lab['r_edge'] + 0.015) * math.cos(a), (lab['r_edge'] + 0.015) * math.sin(a)),
                         ((r - 0.02) * math.cos(a), (r - 0.02) * math.sin(a))]


# -- scoring ------------------------------------------------------------------------------------

def _in_ring(x, y, L, margin=0.0):
    """True if (x, y) is on a ring, or within margin of one."""
    r = math.hypot(x, y)
    th = math.degrees(math.atan2(y, x))
    for ring, a, b, _ in L['wedges']:
        if RADII[ring] - 0.01 - margin <= r <= RADII[ring] + RING_WIDTH + 0.01 + margin:
            pad = math.degrees(margin / max(r, 1e-6))
            if in_span(th, min(a, b) - pad, max(a, b) + pad):
                return True
    return False


def _overlap(b1, b2, pad=0.01):
    return not (b1[2] + pad <= b2[0] or b2[2] + pad <= b1[0] or b1[3] + pad <= b2[1] or b2[3] + pad <= b1[1])


def content_box(L):
    """(x0, y0, x1, y1) of rings, labels and leader lines, in data units (the title is not part of
    the chart's rectangle; render() puts it above)."""
    xs, ys = [], []
    for ring, a, b, _ in L['wedges']:
        R = RADII[ring] + RING_WIDTH
        n = max(2, int(abs(b - a)))
        for i in range(n + 1):
            t = math.radians(a + (b - a) * i / n)
            xs.append(R * math.cos(t))
            ys.append(R * math.sin(t))
        for t in (0, 90, 180, 270):
            if in_span(t, a, b):
                xs.append(R * math.cos(math.radians(t)))
                ys.append(R * math.sin(math.radians(t)))
    for lab in L['outside']:
        x0, y0, x1, y1 = text_box(lab)
        xs += [x0, x1]
        ys += [y0, y1]
        for px, py in lab['leader']:
            xs.append(px)
            ys.append(py)
    return min(xs), min(ys), max(xs), max(ys)


def _segments_cross(p1, p2, p3, p4):
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2 = orient(p3, p4, p1), orient(p3, p4, p2)
    d3, d4 = orient(p1, p2, p3), orient(p1, p2, p4)
    return d1 * d2 < 0 and d3 * d4 < 0


def _segment_hits_box(p, q, b, n=12):
    return any(b[0] <= p[0] + (q[0] - p[0]) * i / n <= b[2] and b[1] <= p[1] + (q[1] - p[1]) * i / n <= b[3]
               for i in range(n + 1))


def score(L, partial=False):
    """(area in square inches, mean tilt of label letters from upright in degrees), or None if a
    label or leader line touches a ring, another label or another leader line."""
    if L.get('invalid'):
        return None
    if partial:                                       # only the labels placed so far
        return score(dict(L, outside=[lab for lab in L['outside'] if 'leader' in lab]))
    boxes = [text_box(lab) for lab in L['outside']]
    legs = [list(zip(lab['leader'], lab['leader'][1:])) for lab in L['outside']]
    for segs in legs:                                 # leader lines must not cross a ring...
        for p, q in segs:
            for i in range(1, 21):
                if _in_ring(p[0] + (q[0] - p[0]) * i / 20, p[1] + (q[1] - p[1]) * i / 20, L):
                    return None
    for (i, si), (j, sj) in itertools.combinations(enumerate(legs), 2):   # ...each other...
        if any(_segments_cross(p, q, r, t) for p, q in si for r, t in sj):
            return None
    for i, segs in enumerate(legs):                   # ...or a label (its own only at the end point)
        for j, b in enumerate(boxes):
            if i != j and any(_segment_hits_box(p, q, b) for p, q in segs):
                return None
        if _crosses_own_label(L['outside'][i]):
            return None
    for b in boxes:                                   # labels must not touch a ring...
        for i in range(5):
            for j in range(3):
                if _in_ring(b[0] + (b[2] - b[0]) * i / 4, b[1] + (b[3] - b[1]) * j / 2, L):
                    return None
    for b1, b2 in itertools.combinations(boxes, 2):   # ...or each other
        if _overlap(b1, b2):
            return None
    x0, y0, x1, y1 = content_box(L)
    area = (x1 - x0) * (y1 - y0) * (PT_PER_UNIT / 72) ** 2
    num = den = 0.0
    for lab in L['arc']:
        for _, w, _, _, r in lab['chars']:
            num += w * tilt(r)
            den += w
    for lab in L['outside']:
        den += sum(text_width_pt(s, LABEL_PT) for s in lab['text'].split('\n')) * UNIT
    return area, num / den


def search_layouts(data, title, starts=range(0, 360, 5), max_tilts=(180, 60, 40, 25, 10), fixed=None):
    """Score every combination of parameters (with the fixed ones, e.g. label_style, held);
    returns [(area, tilt, params)] of valid layouts."""
    out = []
    for start, cw, o1, o2, o3, o4, mt, polar in itertools.product(
            starts, (True, False), itertools.permutations(LEVELS[1][1]), (LEVELS[2][1], ('C-', 'C+', 'D')),
            itertools.permutations(LEVELS[3][1]), (LEVELS[4][1], LEVELS[4][1][::-1]), max_tilts, (False, True)):
        params = {'start': float(start), 'clockwise': cw, 'order1': o1, 'order2': o2, 'order3': o3,
                  'order4': o4, 'max_tilt': float(mt), 'polar_labels': polar, **(fixed or {})}
        s = score(layout(data, title, params))
        if s:
            out.append((s[0], s[1], params))
    return out


# -- drawing ------------------------------------------------------------------------------------

def render(L, data, stem, dpi=600, palette='published'):
    """Draw a layout with its legend to stem.{svg,pdf,png}; returns the paths."""
    colors = PALETTES[palette]
    x0, y0, x1, y1 = content_box(L)
    # title centred over the chart, above it; the legend's top is level with the title's top
    tw = text_width_pt(L['title'], TITLE_PT, True) * UNIT
    tx, ty = (x0 + x1) / 2, y1 + 0.14
    top = ty + CAP * TITLE_PT * UNIT
    legend_w = max(text_width_pt(h, LEGEND_PT) for h, _ in LEGEND) * UNIT + 0.1
    lx = x1 + 0.22
    bottom = min(y0, top - legend_height(data))
    x0 = min(x0, tx - tw / 2)
    W = (lx + legend_w - x0) + 0.1
    H = (top - bottom) + 0.1
    y0 = bottom
    fig = plt.figure(figsize=(W * PT_PER_UNIT / 72, H * PT_PER_UNIT / 72), facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0 - 0.05, x0 - 0.05 + W)
    ax.set_ylim(y0 - 0.05, y0 - 0.05 + H)
    ax.set_aspect('equal')
    ax.axis('off')
    for ring, a, b, k in L['wedges']:
        lo, hi = min(a, b), max(a, b)
        ax.add_patch(Wedge((0, 0), RADII[ring] + RING_WIDTH, lo, hi, width=RING_WIDTH, facecolor=colors[k],
                           edgecolor='black', linewidth=EDGE_PT))
    ax.add_patch(Wedge((0, 0), RADII[1], 0, 360, facecolor='white', edgecolor='none'))
    for lab in L['arc']:
        for ch, _, x, y, rot in lab['chars']:
            ax.text(x, y, ch, fontsize=lab['size'], family=FONT, rotation=rot, rotation_mode='anchor',
                    ha='center', va='baseline')
    for lab in L['outside']:
        xs, ys = zip(*lab['leader'])
        ax.plot(xs, ys, color='black', linewidth=0.5, solid_joinstyle='round', solid_capstyle='round')
        b = text_box(lab)                                # vertical centre after any attachment shift
        ax.text(lab['x'], (b[1] + b[3]) / 2, lab['text'], fontsize=LABEL_PT, family=FONT, ha=lab['ha'], va='center',
                linespacing=1.15, multialignment='center' if lab['ha'] == 'center' else lab['ha'])
    ax.text(0, 0, L['center'], fontsize=LEGEND_PT + 1, family=FONT, ha='center', va='center')
    ax.text(tx, ty, L['title'], fontsize=TITLE_PT, family=FONT, fontweight='bold', ha='center', va='baseline')
    legend(ax, data, lx, top, colors)
    paths = []
    for ext in ('svg', 'pdf', 'png'):
        path = Path(f'{stem}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white', metadata={'Date': None} if ext == 'svg' else None)
        paths.append(path)
    plt.close(fig)
    return paths


def legend_blocks(data):
    present = {k for lvl in ('level1', 'level2', 'level3', 'level4') for k, v in data[lvl].items() if v}
    blocks = [(h, [(k, s) for k, s in items if k in present]) for h, items in LEGEND]
    return [b for b in blocks if b[1]]


def legend_height(data):
    blocks = legend_blocks(data)
    n_lines = sum(1 + len(items) for _, items in blocks) + 0.6 * (len(blocks) - 1)
    return n_lines * 1.35 * LEGEND_PT * UNIT


def legend(ax, data, x, top, colors=COLORS):
    """The legend with its first header's top at y = top."""
    sw = 7.0 * UNIT
    pitch = 1.35 * LEGEND_PT * UNIT
    blocks = legend_blocks(data)
    cap = CAP * LEGEND_PT * UNIT
    y = top - cap                                     # first header baseline
    for header, items in blocks:
        ax.text(x, y, header, fontsize=LEGEND_PT, family=FONT, ha='left', va='baseline')
        y -= pitch
        for k, s in items:
            ax.add_patch(Rectangle((x, y + cap / 2 - sw / 2), sw, sw, facecolor=colors[k], edgecolor='black',
                                   linewidth=0.4, clip_on=False))
            ax.text(x + sw + 3.5 * UNIT, y, s, fontsize=LEGEND_PT, family=FONT, ha='left', va='baseline')
            y -= pitch
        y -= 0.6 * pitch


# -- outputs ------------------------------------------------------------------------------------

def counts_rows(tag, data):
    rows = []
    parents = {'level1': data['total'], 'level2': sum(data['level2'].values()),
               'level3': sum(data['level3'].values()), 'level4': sum(data['level4'].values())}
    ring = {'level1': 'inner: array count', 'level2': 'second: array type (1 array)',
            'level3': 'third: pLAM (4qA)', 'level4': 'outer: repeat units (functional pLAM)'}
    for lvl in ('level1', 'level2', 'level3', 'level4'):
        for k, v in data[lvl].items():
            n = parents[lvl]
            rows.append((tag, ring[lvl], NAMES.get(k, k), v, n, f'{100 * v / n:.1f}' if n else ''))
    return rows


def datasets(cfg):
    """[(tag, counts, title)] for all haplotypes and for those passing the assembly QC."""
    df = pd.read_csv(cfg['inputs']['analysis_tsv'], sep='\t')
    plam = pd.read_csv(cfg['inputs']['plam_tsv'], sep='\t')
    qc = pd.read_csv(cfg['inputs']['qc_classification_tsv'], sep='\t')
    chrom = cfg.get('sunburst', {}).get('chromosome', 'chr4')
    passing = set(qc.loc[qc['qc_pass'].astype(str) == 'True', 'contig'])
    df_qc = df[df['haplotype'].str.split('_', n=2).str[2].isin(passing)]
    out = []
    for tag, frame in (('all', df), ('qc', df_qc)):   # the QC filter is described in the caption
        data = pad.calculate_distributions(frame, chrom, plam_df=plam)
        out.append((tag, data, f'{chrom} (n={data["total"]})'))
    return chrom, out


def params_from_cfg(cfg):
    p = dict(cfg.get('sunburst', {}).get('layout', {}))
    for k in ('order1', 'order2', 'order3', 'order4', 'inside'):
        if k in p:
            p[k] = tuple(p[k])
    return p


def build(cfg, out):
    """Write both versions of the chr4 sunburst and the counts table; returns the paths."""
    out = Path(out)
    chrom, sets = datasets(cfg)
    params = params_from_cfg(cfg)
    paths, rows = [], []
    for tag, data, title in sets:
        L = layout(data, title, params)
        if score(L) is None:
            raise ValueError(f'sunburst ({tag}): a label collides with a ring or another label; '
                             'choose other layout parameters (sunburst: layout: in the config)')
        stem = out / (f'sunburst_{chrom}' + ('' if tag == 'all' else '_qc'))
        paths += render(L, data, stem, palette=cfg.get('sunburst', {}).get('palette', 'published'))
        rows += counts_rows(tag, data)
    tsv = out / f'sunburst_{chrom}_counts.tsv'
    pd.DataFrame(rows, columns=['version', 'ring', 'category', 'n', 'of', 'percent']).to_csv(tsv, sep='\t', index=False)
    return paths + [tsv]
