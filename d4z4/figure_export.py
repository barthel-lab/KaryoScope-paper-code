#!/usr/bin/env python3
"""
Export track views of chosen haplotypes in the Supplementary Figure S10 style, plus colour
legends, as separate files for use as panels in other figures (see make_figures.py --export).

Every file is drawn at the S10 scale (track width, fonts, line widths), so track views and
legends keep matching text sizes when placed at 100%. By default all exported haplotypes share
one kb scale (the widest haplotype's window), so array lengths can be compared across files;
--scale own gives each its own window, as in S10. Crowded element numbers are thinned with the
rule applied by hand in the published figures (keep the array's total, drop the label next to it),
and the label interval is widened where a shared scale draws an array compressed.

Track views, for haplotypes given as contig IDs (sample#hap#accession):
    <sample>_hap<N>                  one haplotype, with row labels
    <name>_side_by_side              all haplotypes in one row; row labels on the first only (as in S10)
    <name>_stacked                   one haplotype per row, each with row labels
Colour legends, listing only the features drawn in the exported tracks (or all of them):
    legend_wide                      one column per legend group, groups of more than 6 entries
                                     split in two (as in S10)
    legend_<k>col                    the groups packed into k columns
Each file is written as .svg (text kept as text), .pdf and .png.

Legend packing: the groups are stacked into k columns, each column as wide as its widest
group. A group with many entries may be split into side-by-side sub-columns (Repeats 5 + 5 in
S10), but only into sub-columns that are all full (10 entries as 5 + 5, never 5 as 3 + 2), so
no split leaves a hole. Of all ways to assign the groups to k columns (keeping the S10 order
within each column) and to split them, the one with the smallest bounding box is drawn.
"""

import itertools
import re
from pathlib import Path

import figure_build as fb
import figure_compose as fc
import figure_data as fd
import plot_d4z4_arrays as pda

# Page geometry, in points
MARGIN = 4.0          # white space around every page
ROW_GAP = 8.0         # between stacked haplotypes
# Legend geometry, in points (the item sizes and pitch come from figures_layout.yaml track_legend)
HEADER_TO_ITEMS = 16.35  # header baseline to first item baseline (S10)
HEADER_DX = -7.4         # header x relative to the item text x (S10)
SUBCOL_GAP = 5.0         # end of a label to the next sub-column's swatch, within a group
COL_GAP = 8.0            # between legend columns
BLOCK_GAP = 8.0          # between groups stacked in one column
CAP, DESC = 0.716, 0.212  # Arial cap height and descent, as fractions of the font size


def label(contig):
    """'HG00621#1#CM087849.1' -> 'HG00621_hap1'."""
    sample, hap, _ = contig.split('#')
    return f'{sample}_hap{hap}'


def group_name(contigs):
    samples = {c.split('#')[0] for c in contigs}
    return samples.pop() if len(samples) == 1 else '_'.join(label(c) for c in contigs)


# -- track views ------------------------------------------------------------------------------

def _page(name, rows, cfg, layout, out, dpi, window_size, width_pt=None):
    """Draw rows of haplotypes (a list of lists of contigs) on one page in the S10 geometry.

    width_pt: the page width; the tracks are lengthened (or shortened) to fill it while text and
    line widths keep their S10 point sizes (the drawing canvas widens, its scale does not change)."""
    S10 = layout['figS10']
    tv = layout['track_view']
    (c0, c1), (d0, _) = S10['track_columns'][0], S10['track_columns'][1]
    col_w, col_pitch = c1 - c0, d0 - c0
    if width_pt and len(rows[0]) == 1 and max(len(r) for r in rows) == 1:
        s0 = col_w / (tv['width'] - tv['left_margin'] - tv['right_margin'])
        new_w = width_pt - (c0 + tv['right_margin'] * s0 + 2 * MARGIN)
        tv = dict(tv, width=tv['left_margin'] + tv['right_margin'] + new_w / s0)
        layout = dict(layout, track_view=tv)
        col_w = new_w
    ncols = max(len(r) for r in rows)
    columns = {i: [c0 + i * col_pitch, c0 + col_w + i * col_pitch] for i in range(ncols)}
    tracks = [{'panel': f'R{r}', 'column': i, 'ruler_y': 0.0} for r, row in enumerate(rows) for i in range(len(row))]
    layout[name] = dict(S10, letters={}, track_legend=[], track_columns=columns, tracks=tracks)
    cfg[name] = {f'R{r}': {'haplotypes': list(row)} for r, row in enumerate(rows)}
    ctx = fb.FigureContext(name, cfg, layout, out)
    ctx.track_style = dict(ctx.track_style, auto_hide_crowded_labels=True, label_min_gap_em=2.4,
                           window_size=window_size)
    s = col_w / (tv['width'] - tv['left_margin'] - tv['right_margin'])
    title_y = ctx.track_style.get('title_font') or max(9, int(int(15 * tv['height_ratio']) * 0.85))
    axis_off = (title_y + S10['title_axis_gap']) * s   # ruler distance below the drawing's top
    y = bottom = MARGIN
    tops = []
    for r, row in enumerate(rows):
        if r:
            y = bottom + ROW_GAP
        tops.append(y)
        for i, contig in enumerate(row):
            ctx.track_view(contig, i, y + axis_off, f'R{r}_{i}')
            h = float(re.search(r'\bheight="([\d.]+)"', ctx.items[-1]).group(1))
            bottom = max(bottom, y + h)
    # In S10 page coordinates the row labels use x = 0 to c0; the last track's drawing ends
    # right_margin drawing units after its track box
    right = c0 + (ncols - 1) * col_pitch + col_w + tv['right_margin'] * s
    ctx.items = [f'<g transform="translate({MARGIN:.3f},0)">'] + ctx.items + ['</g>']
    w, h = round(right + 2 * MARGIN, 2), round(bottom + MARGIN, 2)
    ctx.L['size_pt'] = [w, h]
    ctx.L['png_px'] = [round(w / 72 * dpi), round(h / 72 * dpi)]
    paths = ctx.write()
    del layout[name], cfg[name]
    return paths, ctx.renderer.used_colors, tops


def shared_window(cfg, layout, contigs):
    """The largest adaptive window among the haplotypes, in bp."""
    from figure_tracks import TrackRenderer
    df = fd.with_config(cfg)
    renderer = TrackRenderer(cfg, df, {})
    return max(renderer.window_size(fd.haplotype_id(df, c)) for c in contigs)


def export_tracks(cfg, layout, contigs, out, arrangements, dpi, scale='shared'):
    """Write the requested track-view arrangements; returns {file stem: paths} and the colours used."""
    window_size = shared_window(cfg, layout, contigs) if scale == 'shared' else None
    pages = []
    if 'single' in arrangements:
        pages += [(label(c), [[c]]) for c in contigs]
    if len(contigs) > 1 and 'side_by_side' in arrangements:
        pages.append((f'{group_name(contigs)}_side_by_side', [list(contigs)]))
    if len(contigs) > 1 and 'stacked' in arrangements:
        pages.append((f'{group_name(contigs)}_stacked', [[c] for c in contigs]))
    written, used = {}, {}
    for name, rows in pages:
        paths, u, _ = _page(name, rows, cfg, layout, out, dpi, window_size)
        written[name] = paths
        for kind, feats in u.items():
            used.setdefault(kind, set()).update(feats)
    return written, used


# -- legends ----------------------------------------------------------------------------------

def legend_groups(layout, used=None, drawn_fills=None):
    """[(header, [(label, key), ...]), ...] in S10 order. Header-less sections continue the
    previous group (S10 splits Repeats over two sections). With used/drawn_fills, entries for
    features not drawn are left out, and so are groups left empty."""
    aliases = {'arm_multigroup1': {'arm_multigroup1', 'p_arm_specific', 'q_arm_specific'}}
    groups = []
    for sec in layout['figS10']['track_legend']:
        items = []
        for lbl, key in sec['items']:
            kind, _, feat = key.partition(':')
            if used is not None:
                if kind in ('gene', 'region', 'repeat', 'subtelomeric') and \
                        not aliases.get(feat, {feat}) & set(used.get(kind, ())):
                    continue
                if kind == 'novel' and not any('novel' in v for v in used.values()):
                    continue
                if kind in ('d4z4', 'triangle') and getattr(pda, feat).lower() not in drawn_fills:
                    continue
            items.append((lbl, key))
        if sec['header'] is None and groups:
            groups[-1][1].extend(items)
        elif items:
            groups.append((sec['header'], items))
    return [g for g in groups if g[1]]


def set_partitions(n, k):
    """Every way to put items 0..n-1 into k non-empty columns, each column in increasing order,
    columns ordered by their first item (restricted growth strings)."""
    def grow(i, labels, used):
        if i == n:
            if used == k:
                yield [[j for j in range(n) if labels[j] == c] for c in range(k)]
            return
        if n - i < k - used:     # not enough items left to fill the remaining columns
            return
        for c in range(min(used + 1, k)):
            labels.append(c)
            yield from grow(i + 1, labels, max(used, c + 1))
            labels.pop()
    yield from grow(0, [], 0)


class LegendPacker:
    def __init__(self, layout, measure):
        self.spec = layout['track_legend']
        self.measure = measure   # (text, size, bold) -> width in points

    def block(self, header, items, nsub):
        """Size of one group drawn with its entries in nsub sub-columns, plus the sub-column layout."""
        sp = self.spec
        per = -(-len(items) // nsub)
        cols = [items[i:i + per] for i in range(0, len(items), per)]
        text_x = -HEADER_DX                                   # block's left edge is the header's left edge
        xs, x = [], text_x
        for c in cols:
            xs.append(x)
            x += max(self.measure(lbl, sp['item_size'], False) for lbl, _ in c) + SUBCOL_GAP - sp['swatch_dx']
        items_right = x - SUBCOL_GAP + sp['swatch_dx']
        header_right = self.measure(header or '', sp['header_size'], True)
        hb = CAP * sp['header_size']                          # header baseline below the block top
        rows = max(len(c) for c in cols)
        height = hb + HEADER_TO_ITEMS + (rows - 1) * sp['item_pitch'] + DESC * sp['item_size']
        return {'header': header, 'cols': cols, 'xs': xs, 'hb': hb,
                'w': max(items_right, header_right), 'h': height, 'nsub': len(cols)}

    def pack(self, groups, k, wide=False, max_width=None, max_height=None, minimize='area'):
        """Best packing of the groups into k columns: (width, height, [(block, x, y), ...]), or None
        if no packing fits max_width/max_height.

        wide: one column per group, and only groups of more than 6 entries split, in two (S10).
        minimize: 'area', 'height' (for a legend below a figure) or 'width' (beside it)."""
        options = []
        for header, items in groups:
            if wide:
                options.append([self.block(header, items, 1 if len(items) <= 6 else 2)])
                continue
            # only splits into equally full sub-columns: no holes
            options.append([self.block(header, items, n) for n in range(1, min(3, len(items)) + 1)
                            if len(items) % n == 0])
        n = len(groups)
        best = None
        for spans in set_partitions(n, k):
            # reading order (down each column, then across) against the S10 order, for tie-breaks
            order = [j for sp in spans for j in sp]
            inversions = sum(a > b for i, a in enumerate(order) for b in order[i + 1:])
            for choice in itertools.product(*options):
                width = sum(max(choice[j]['w'] for j in sp) for sp in spans) + COL_GAP * (k - 1)
                height = max(sum(choice[j]['h'] for j in sp) + BLOCK_GAP * (len(sp) - 1) for sp in spans)
                if (max_width and width > max_width) or (max_height and height > max_height):
                    continue
                splits = sum(b['nsub'] - 1 for b in choice)
                # areas within 50 pt^2 (well under 1%) count as equal: then fewer splits, then
                # closer to the S10 order
                primary = {'area': round(width * height / 50), 'height': round(height, 1),
                           'width': round(width, 1)}[minimize]
                key = (primary, round(width * height / 50), splits, inversions)
                if best is None or key < best[0]:
                    best = (key, width, height, spans, choice)
        if best is None:
            return None
        _, width, height, spans, choice = best
        placed, x = [], 0.0
        for sp in spans:
            y = 0.0
            for j in sp:
                placed.append((choice[j], x, y))
                y += choice[j]['h'] + BLOCK_GAP
            x += max(choice[j]['w'] for j in sp) + COL_GAP
        return width, height, placed

    def svg_items(self, placed, ox, oy, sans, color_maps):
        sp = self.spec
        out = []
        for b, bx, by in placed:
            x0, y0 = ox + bx, oy + by
            if b['header']:
                out.append(fb._txt(b['header'], x0, y0 + b['hb'], sp['header_size'], sans, bold=True))
            for col, cx in zip(b['cols'], b['xs']):
                for r, (lbl, key) in enumerate(col):
                    tx, ty = x0 + cx, y0 + b['hb'] + HEADER_TO_ITEMS + r * sp['item_pitch']
                    sx, sy = tx + sp['swatch_dx'], ty + sp['swatch_dy']
                    kind, _, feat = key.partition(':')
                    a = sp['swatch']
                    if kind == 'triangle':
                        out.append(f'<path d="M{sx + a / 2:.3f},{sy:.3f} L{sx:.3f},{sy + a:.3f} '
                                   f'L{sx + a:.3f},{sy + a:.3f} Z" fill="{getattr(pda, feat)}"/>')
                    else:
                        color = (pda.DEFAULT_COLOR if kind == 'novel' else
                                 getattr(pda, feat) if kind == 'd4z4' else color_maps[kind][feat])
                        out.append(fb._rect(sx, sy, a, a, color, '#000000', 0.5))
                    out.append(fb._txt(lbl, tx, ty, sp['item_size'], sans))
        return out


def text_measure():
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath
    t2p = TextToPath()

    def measure(s, size, bold):
        if not s:
            return 0.0
        fp = FontProperties(family=['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'], size=size,
                            weight='bold' if bold else 'normal')
        return t2p.get_text_width_height_descent(s, fp, ismath=False)[0]
    return measure


def export_legends(cfg, layout, out, columns, dpi, used=None, drawn_fills=None):
    """Write one legend per entry of columns ('wide' or a number of columns)."""
    groups = legend_groups(layout, used, drawn_fills)
    packer = LegendPacker(layout, text_measure())
    color_maps = pda.load_all_colors(cfg['inputs']['colors_dir'])
    sans = layout['fonts']['sans']
    written, seen = {}, {}
    for c in columns:
        k = len(groups) if c == 'wide' else min(int(c), len(groups))
        name = 'legend_wide' if c == 'wide' else f'legend_{k}col'
        if name in seen:
            continue
        seen[name] = k
        width, height, placed = packer.pack(groups, k, wide=(c == 'wide'))
        items = packer.svg_items(placed, MARGIN, MARGIN, sans, color_maps)
        W, H = round(width + 2 * MARGIN, 2), round(height + 2 * MARGIN, 2)
        stem = Path(out) / name
        fc.compose(items, W, H, f'{stem}.svg')
        fc.render(f'{stem}.svg', f'{stem}.pdf', f'{stem}.png', png_size=[round(W / 72 * dpi), round(H / 72 * dpi)])
        ink = sum(b['w'] * b['h'] for b, _, _ in placed)
        written[name] = (W, H, 1 - ink / (width * height))
    return written


def export(cfg, layout, contigs, out, arrangements, legend_columns, dpi, legend_features, scale='shared'):
    out = Path(out)
    (out / 'panels').mkdir(parents=True, exist_ok=True)
    df = fd.load_analysis(cfg)
    for c in contigs:
        fd.haplotype_id(df, c)   # every contig must be in the analysis table (raises otherwise)
    tracks, used = export_tracks(cfg, layout, contigs, out, arrangements, dpi, scale)
    drawn_fills = {m.lower() for paths in tracks.values() for p in paths if p.suffix == '.svg'
                   for m in re.findall(r'fill="(#[0-9A-Fa-f]{6})"', p.read_text())}
    if legend_features == 'all':
        used = drawn_fills = None
    legends = export_legends(cfg, layout, out, legend_columns, dpi, used, drawn_fills)
    return tracks, legends


# -- composed figures (HPRC2 Fig. 6f, Supplementary Fig. 42) ----------------------------------

def _svg_size(path):
    head = Path(path).read_text()[:400]
    w, h = (float(re.search(rf'\b{k}="([\d.]+)pt"', head).group(1)) for k in ('width', 'height'))
    return w, h


def composed_figure(name, spec, cfg, layout, out, dpi=600):
    """One figure: the haplotypes stacked on one kb scale, optional panel letters, and the colour
    legend packed beside ('right') or below ('below') the tracks so that it fits their height or
    width. spec is a figures_config.yaml entry: haplotypes, letters, legend, legend_features, and
    width_pt (page width; the tracks lengthen to fill it, text keeps its size)."""
    import figure_compose as fc
    out = Path(out)
    work = out / 'panels'
    (work / 'panels').mkdir(parents=True, exist_ok=True)
    contigs = spec['haplotypes']
    df = fd.load_analysis(cfg)
    for c in contigs:
        fd.haplotype_id(df, c)   # raises for a contig that is not in the analysis table
    window = shared_window(cfg, layout, contigs) if spec.get('scale', 'shared') == 'shared' else None
    paths, used, tops = _page(f'{name}_tracks', [[c] for c in contigs], cfg, layout, work, dpi, window,
                              width_pt=spec.get('width_pt'))
    tracks_svg = paths[0]
    W, H = _svg_size(tracks_svg)
    drawn = {m.lower() for m in re.findall(r'fill="(#[0-9A-Fa-f]{6})"', tracks_svg.read_text())}
    if spec.get('legend_features', 'drawn') == 'all':
        used = drawn = None
    groups = legend_groups(layout, used, drawn)
    packer = LegendPacker(layout, text_measure())
    where = spec.get('legend', 'below')
    ks = range(1, len(groups) + 1)
    if where == 'right':        # as narrow as possible while no taller than the tracks
        fit = [o for o in (packer.pack(groups, k, max_height=H - 2 * MARGIN, minimize='width') for k in ks) if o]
        lw, lh, placed = min(fit, key=lambda o: (o[0], o[1])) if fit else \
            min((packer.pack(groups, k, minimize='height') for k in ks), key=lambda o: o[1])
        ox, oy = W + 6.0, MARGIN
    else:                       # as short as possible while no wider than the tracks
        fit = [o for o in (packer.pack(groups, k, max_width=W - 2 * MARGIN - 8.0, minimize='height')
                           for k in ks) if o]
        lw, lh, placed = min(fit, key=lambda o: (o[1], o[0])) if fit else \
            min((packer.pack(groups, k, minimize='width') for k in ks), key=lambda o: o[0])
        ox, oy = MARGIN + 8.0, H + 2.0
    color_maps = pda.load_all_colors(cfg['inputs']['colors_dir'])
    items = [fc.nest_svg(tracks_svg, 0, 0, W, H, prefix='T_')]
    items += packer.svg_items(placed, ox, oy, layout['fonts']['sans'], color_maps)
    size = spec.get('letter_size', 10)
    for letter, top in zip(spec.get('letters', []), tops):
        items.append(fb._txt(letter, MARGIN, top + CAP * size, size, layout['fonts']['sans'], bold=True))
    FW = round(max(W, ox + lw + MARGIN), 2)
    FH = round(max(H, oy + lh + MARGIN), 2)
    stem = out / name
    fc.compose(items, FW, FH, f'{stem}.svg')
    fc.render(f'{stem}.svg', f'{stem}.pdf', f'{stem}.png', png_size=[round(FW / 72 * dpi), round(FH / 72 * dpi)])
    return [Path(f'{stem}.{ext}') for ext in ('svg', 'pdf', 'png')]
