#!/usr/bin/env python3
"""
Build Figure 4, Supplementary Figure S9 and Supplementary Figure S10 (see make_figures.py).

Track views are drawn by plot_d4z4_arrays.py and placed so their track boxes span the
published columns; the scatter plots are matplotlib; the waffles, diplotype tables, colour
legends, schematic and panel letters are drawn directly as SVG at their published positions.
"""

import copy
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex

import figure_compose as fc
import figure_data as fd
import figure_schematic
import plot_d4z4_arrays as pda
from create_d4z4_figure import C_BLUE, C_EMERALD, C_GRAY, CONFIG_ALPHA, CONFIG_COLORS
from figure_tracks import TrackRenderer

PT = 1 / 72


def _txt(s, x, y, size, font, bold=False, anchor='start', fill='#000000', rotate=None):
    return fc.text_el(s, x, y, size, font, 'bold' if bold else 'normal', anchor, fill, rotate)


def _rect(x, y, w, h, fill, stroke=None, sw=0, opacity=None):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ''
    o = f' fill-opacity="{opacity}"' if opacity is not None and opacity < 1 else ''
    return f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" fill="{fill}"{o}{s}/>'


def omit_panels(L, omitted):
    """Layout without the omitted panels (a stack of full-width panels, as in S10).

    Each omitted panel's tracks and letter are dropped, everything below it moves up by the
    panel's height (its letter to the next letter), and the page gets shorter by the same
    amount. Panels keep their slot names in the config; the printed letters run on from A.
    """
    if not omitted:
        return L
    L = copy.deepcopy(L)
    order = sorted(L['letters'], key=lambda k: L['letters'][k][1])
    for p in sorted(omitted, key=lambda k: -L['letters'][k][1]):    # bottom-up, so y's stay valid
        y0 = L['letters'][p][1]
        later = [k for k in order if L['letters'][k][1] > y0 and k not in omitted]
        if not later:
            raise ValueError(f'omit_panels: {p} is the last panel; only inner panels can be omitted')
        dy = L['letters'][later[0]][1] - y0
        L['tracks'] = [t for t in L['tracks'] if t['panel'] != p]
        for t in L['tracks']:
            if t['ruler_y'] > y0:
                t['ruler_y'] -= dy
        for k in later:
            L['letters'][k][1] -= dy
        for sec in L.get('track_legend', []):
            for key in ('header_at', 'items_at'):
                if sec[key] and sec[key][1] > y0:
                    sec[key][1] -= dy
        w, h = L['size_pt']
        L['png_px'] = [L['png_px'][0], round(L['png_px'][1] * (h - dy) / h)]
        L['size_pt'] = [w, h - dy]
        del L['letters'][p]
    kept = [k for k in order if k not in omitted]
    L['printed_letters'] = {k: chr(ord('a') + i) for i, k in enumerate(kept)}
    return L


class FigureContext:
    def __init__(self, name, cfg, layout, out):
        self.name, self.cfg, self.layout = name, cfg, layout
        self.L = omit_panels(layout[name], cfg[name].get('omit_panels', []))
        self.out = Path(out)
        self.panels = self.out / 'panels'
        self.items = []
        self.sans = layout['fonts']['sans']
        self.plot_font = layout['fonts']['plot']
        self.df = fd.with_config(cfg)
        tv = layout['track_view']
        self.track_style = dict(tv['style'], font_family=self.sans, **self.L.get('track_style', {}))
        self.renderer = TrackRenderer(cfg, self.df, self.track_style)
        # output.width_mm: the page is written at that width (the layout stays in published-page
        # points), so every size in the layout prints at `scale` times its value
        w_pt = self.L['size_pt'][0]
        width_mm = layout.get('output', {}).get('width_mm')
        self.scale = width_mm / (w_pt / 72 * 25.4) if width_mm else 1.0

    # -- track views --------------------------------------------------------------------
    def track_view(self, contig, column, ruler_y, tag):
        tv = self.layout['track_view']
        hap = fd.haplotype_id(self.df, contig)
        hidden = self.layout.get('hidden_element_labels', {}).get(contig, [])
        # The drawing is scaled by s so that its track box spans the column
        x0, x1 = self.L['track_columns'][column]
        plot_width = tv['width'] - tv['left_margin'] - tv['right_margin']
        s = (x1 - x0) / plot_width
        # Row labels ('Satellites', ...) are drawn for the left column only; lines (track borders,
        # ruler, ticks, array bars) keep their printed width in points
        self.renderer.style = dict(self.track_style, hidden_element_labels={tuple(h) for h in hidden},
                                   track_labels=(column == 0), title_axis_gap=self.L['title_axis_gap'],
                                   stroke_width=tv['line_width_pt'] / s)
        path = self.panels / f'{self.name}_{tag}.svg'
        w, h = self.renderer.draw(hap, str(path), tv['width'], tv['left_margin'], tv['right_margin'],
                                  tv['height_ratio'])
        # the drawing's ruler sits at (title baseline + title_axis_gap) drawing units from its top
        title_y = self.renderer.style.get('title_font') or max(9, int(int(15 * tv['height_ratio']) * 0.85))
        axis_y = title_y + self.L['title_axis_gap']
        self.items.append(fc.nest_svg(path, x0 - tv['left_margin'] * s, ruler_y - axis_y * s,
                                      w * s, h * s, prefix=f'{tag}_'))

    def all_tracks(self):
        used = {}
        for t in self.L['tracks']:
            k = used.get(t['panel'], 0)
            used[t['panel']] = k + 1
            listed = self.cfg[self.name][t['panel']]['haplotypes']
            if k >= len(listed):    # fewer haplotypes than slots: the slot stays empty
                continue
            contig = listed[k]
            self.track_view(contig, t['column'], t['ruler_y'], f"{t['panel']}{k + 1}")

    # -- colour legend shared by the track views ----------------------------------------
    def track_legend(self):
        spec = self.layout['track_legend']
        cmaps = self.renderer.color_maps
        covered = set()
        for sec in self.L['track_legend']:
            if sec['header']:
                hx, hy = sec['header_at']
                self.items.append(_txt(sec['header'], hx, hy, spec['header_size'], self.sans, bold=True))
            ix, iy = sec['items_at']
            for k, (label, key) in enumerate(sec['items']):
                y = iy + k * spec['item_pitch']
                sx, sy = ix + spec['swatch_dx'], y + spec['swatch_dy']
                kind, _, name = key.partition(':')
                if kind == 'triangle':
                    color = getattr(pda, name)
                    a = spec['swatch']
                    self.items.append(f'<path d="M{sx + a / 2:.3f},{sy:.3f} L{sx:.3f},{sy + a:.3f} '
                                      f'L{sx + a:.3f},{sy + a:.3f} Z" fill="{color}"/>')
                else:
                    if kind == 'novel':
                        color = pda.DEFAULT_COLOR
                    elif kind == 'd4z4':
                        color = getattr(pda, name)
                    else:
                        color = cmaps[kind][name]
                        covered.add((kind, name))
                    self.items.append(_rect(sx, sy, spec['swatch'], spec['swatch'], color, '#000000', 0.5))
                self.items.append(_txt(label, ix, y, spec['item_size'], self.sans))
        # Every colour drawn in a track must be explained by the legend
        aliases = {('region', 'p_arm_specific'): ('region', 'arm_multigroup1'),
                   ('region', 'q_arm_specific'): ('region', 'arm_multigroup1')}
        for kind, feats in self.renderer.used_colors.items():
            for feat in feats:
                if feat == 'novel':
                    continue
                key = aliases.get((kind, feat), (kind, feat))
                if key not in covered:
                    warnings.warn(f'{self.name}: legend has no entry for {kind} feature {feat}')

    def letters(self):
        fonts = self.layout['fonts']
        size = fonts['letter_size']
        if 'output' in self.layout and fonts.get('letter_size_final'):
            size = round(fonts['letter_size_final'] / self.scale, 3)
        printed = self.L.get('printed_letters', {})
        for letter, (x, y) in self.L['letters'].items():
            self.items.append(_txt(printed.get(letter, letter.lower()), x, y, size, self.sans, bold=True))

    def write(self):
        w, h = self.L['size_pt']
        stem = self.out / self.name
        png_px = self.L['png_px']
        ppi = self.layout.get('output', {}).get('ppi')
        if ppi:
            png_px = [round(w * self.scale / 72 * ppi), round(h * self.scale / 72 * ppi)]
        fc.compose(self.items, w, h, f'{stem}.svg', scale=self.scale)
        fc.render(f'{stem}.svg', f'{stem}.pdf', f'{stem}.png', png_size=png_px)
        return [Path(f'{stem}.{ext}') for ext in ('svg', 'pdf', 'png')]


# -- Figure 4 -----------------------------------------------------------------------------

def waffles(ctx):
    """Fig. 4B: one cell per haplotype, in configured category order, row by row from the top."""
    B = ctx.L['B']
    counts = fd.waffle_counts(ctx.cfg, ctx.df)
    for chrom, w in counts.items():
        ncols = int(np.ceil(np.sqrt(w['n'] * 1.3)))
        (gx, gy), (px, py) = B['grids'][chrom]['origin'], B['grids'][chrom]['pitch']
        cells = [c for c, n in w['counts'] for _ in range(n)]
        for i, cat in enumerate(cells):
            r, c = divmod(i, ncols)
            ctx.items.append(_rect(gx + c * px, gy + r * py, px * B['cell_fraction'], py * B['cell_fraction'],
                                   CONFIG_COLORS[cat], B['cell_edge'], B['cell_edge_width'],
                                   opacity=CONFIG_ALPHA.get(cat, 1.0)))
        tx, ty = B['titles'][chrom]
        ctx.items.append(_txt(f'{chrom} (n={w["n"]})', tx, ty, B['title_size'], ctx.sans, anchor='middle'))
    # Legend: the configured categories (chr4 list, then chr10 additions) present on either
    # chromosome, as in generate_panel_b_waffle(); filled column by column
    total = {}
    for w in counts.values():
        for c, n in w['counts']:
            total[c] = total.get(c, 0) + n
    cats = []
    for cat_list in ctx.cfg['fig4']['B']['categories'].values():
        cats += [c for c in cat_list if c not in cats and total.get(c, 0) > 0]
    leg = B['legend']
    nrows = len(leg['rows_baseline'])
    for k, cat in enumerate(cats):
        x, y = leg['columns_x'][k // nrows], leg['rows_baseline'][k % nrows]
        sw, sh = leg['swatch']
        ctx.items.append(_rect(x + leg['swatch_dx'], y - sh, sw, sh, CONFIG_COLORS[cat],
                               opacity=CONFIG_ALPHA.get(cat, 1.0)))
        ctx.items.append(_txt(cat, x, y, leg['size'], ctx.sans))


def diplotype_tables(ctx):
    """Fig. 4G: lower-triangle pairing tables (one unordered pair per individual)."""
    G = ctx.L['G']
    cmap = LinearSegmentedColormap.from_list('custom_blues', ['#F0F4FA', '#A8C4E0', '#4A7FB5', '#1A3A5C'])
    tables = fd.diplotype_tables(ctx.cfg, ctx.df)
    f, ts = ctx.sans, G['text_size']
    for chrom, t in tables.items():
        x0, y0, W, H = G['tables'][chrom]['axes']
        cats, m = t['categories'], t['matrix']
        n = len(cats)
        cw, ch = W / n, H / n
        norm = Normalize(vmin=np.nanmin(m), vmax=np.nanmax(m))
        vmax = np.nanmax(m)
        for i in range(n):
            for j in range(i + 1):
                v = m[i, j]
                ctx.items.append(_rect(x0 + j * cw, y0 + i * ch, cw + 0.05, ch + 0.05, to_hex(cmap(norm(v)))))
                if v > 0:
                    ctx.items.append(_txt(str(int(v)), x0 + (j + 0.5) * cw, y0 + (i + 0.5) * ch + 0.36 * G['cell_size'],
                                          G['cell_size'], ctx.plot_font, anchor='middle',
                                          fill='#FFFFFF' if v > vmax * 0.5 else '#000000'))
        ctx.items.append(f'<rect x="{x0:.3f}" y="{y0:.3f}" width="{W:.3f}" height="{H:.3f}" fill="none" '
                         f'stroke="#000000" stroke-width="{G["frame_width"]}"/>')
        tl = G['tick_length']
        for k, cat in enumerate(cats):
            yc = y0 + (k + 0.5) * ch
            xc = x0 + (k + 0.5) * cw
            ctx.items.append(fc.line_el(x0 - tl, yc, x0, yc, width=G['frame_width']))
            ctx.items.append(fc.line_el(xc, y0 + H, xc, y0 + H + tl, width=G['frame_width']))
            ctx.items.append(_txt(cat, x0 - G['ytick_gap'], yc + G['ytick_dy'], ts, f, anchor='end'))
            ax_, ay_ = xc + G['xtick_dx'], y0 + H + G['xtick_dy']
            ctx.items.append(_txt(cat, ax_, ay_, ts, f, anchor='end', rotate=-45))
        tx, ty = G['tables'][chrom]['title']
        ctx.items.append(_txt(f'{chrom} (n={t["n"]} individuals)', tx, ty, ts, f, anchor='middle'))
        ctx.items.append(_txt('Haplotype 1', x0 + W / 2, y0 + H + G['xlabel_gap'][chrom], ts, f, bold=True,
                              anchor='middle'))
        yx, yy = x0 - G['ylabel_gap'][chrom], y0 + H / 2
        ctx.items.append(_txt('Haplotype 2', yx, yy, ts, f, bold=True, anchor='middle', rotate=-90))


def _mpl_axes(box, axes):
    bx, by, bw, bh = box
    ax0, ay0, aw, ah = axes
    fig = plt.figure(figsize=(bw * PT, bh * PT))
    ax = fig.add_axes([(ax0 - bx) / bw, 1 - (ay0 - by + ah) / bh, aw / bw, ah / bh])
    return fig, ax


def _to_axes_frac(box, axes, x, y):
    """Figure points -> axes fraction."""
    ax0, ay0, aw, ah = axes
    return (x - ax0) / aw, 1 - (y - ay0) / ah


def asymmetry(ctx):
    """Fig. 4H: longer vs shorter haplotype repeat count, single-array individuals."""
    H = ctx.L['H']
    spec = ctx.cfg['fig4']['H']
    pairs, axis_max = fd.asymmetry_pairs(ctx.cfg, ctx.df)
    fig, ax = _mpl_axes(H['box'], H['axes'])
    style = {'chr4': (C_BLUE, 'o'), 'chr10': (C_EMERALD, 's')}
    for chrom in spec['chromosomes']:
        sub = pairs[pairs['chromosome'] == chrom]
        color, marker = style[chrom]
        ax.scatter(sub['shorter'], sub['longer'], s=H['marker_size'], alpha=0.5, color=color, marker=marker,
                   edgecolor='white', linewidth=0.3, zorder=3, label=f'{chrom} (n={len(sub)})')
    ax.plot([0, axis_max], [0, axis_max], 'k-', alpha=0.2, linewidth=1)
    x = np.linspace(1, axis_max, 100)
    for fold in spec['fold_lines']:
        ax.plot(x, x * fold, '--', color=C_GRAY, alpha=0.3, linewidth=0.8)
        # label each guide line near the top of the plot, on the side away from the data
        yl = axis_max * 0.93
        left = fold == max(spec['fold_lines'])
        ax.text(yl / fold + (-1.5 if left else 1.5), yl, f'{fold}\u00d7', fontsize=H['fold_size'],
                color=C_GRAY, ha='right' if left else 'left', va='center')
    top = pairs.nlargest(spec['label_top_n'], 'fold')
    for _, row in top.iterrows():
        placed = H['outlier_labels']
        if ctx.cfg.get('filters', {}).get('assembly_qc') and 'outlier_labels_qc' in H:
            placed = H['outlier_labels_qc']
        lx, ly = placed.get(row['sample'], (None, None))
        if lx is None:  # a new outlier (e.g. after a filter change): put the label beside the point
            xytext, coords = (6, 6), 'offset points'
        else:
            # (lx, ly) is the first line's baseline; va='center' centres the two-line block
            xytext, coords = _to_axes_frac(H['box'], H['axes'], lx, ly + 0.29 * H['label_size']), 'axes fraction'
        ax.annotate(f"{row['sample']}\n({row['fold']:.1f}×)", xy=(row['shorter'], row['longer']),
                    xytext=xytext, textcoords=coords, fontsize=H['label_size'], color=C_GRAY,
                    ha='left', va='center', linespacing=1.2,
                    arrowprops=dict(arrowstyle='-', color=C_GRAY, linewidth=0.25, shrinkA=0, shrinkB=1.5))
    ax.set_xlabel('Shorter haplotype (repeat units)', fontsize=H['axis_size'], labelpad=7.4)
    ax.set_ylabel('Longer haplotype (repeat units)', fontsize=H['axis_size'], labelpad=4.2)
    ax.set_title(f'Haplotype asymmetry (single-array, n={len(pairs)})', fontsize=H['title_size'], pad=3)
    ax.tick_params(labelsize=H['tick_size'], length=3, width=0.6, pad=0.7)
    lx0, ly0, lw, lh = H['legend_box']
    fx, fy = _to_axes_frac(H['box'], H['axes'], lx0, ly0)
    ax.legend(fontsize=H['legend_size'], loc='upper left', bbox_to_anchor=(fx, fy), borderaxespad=0,
              markerscale=0.7, handlelength=1.0, handletextpad=0.3, borderpad=0.35, labelspacing=0.2,
              framealpha=1)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, axis_max)
    path = ctx.panels / 'fig4_H.svg'
    fig.savefig(path, format='svg', transparent=True, metadata={'Date': None})
    plt.close(fig)
    bx, by, bw, bh = H['box']
    ctx.items.append(fc.nest_svg(path, bx, by, bw, bh, prefix='H_'))


def build_fig4(ctx):
    ctx.all_tracks()
    ctx.track_legend()
    waffles(ctx)
    diplotype_tables(ctx)
    asymmetry(ctx)
    ctx.letters()


# -- Supplementary Figure S9 ----------------------------------------------------------------

def hmm_scatter(ctx):
    """Fig. S9B: KaryoScope vs HMM repeat-unit counts, identity line dashed."""
    B = ctx.L['B']
    merged, r = fd.hmm_comparison(ctx.cfg, ctx.df)
    fig, ax = _mpl_axes(B['box'], B['axes'])
    ax.scatter(merged['hmm_count'], merged['karyoscope_count'], s=B['marker_size'], alpha=0.6,
               color=B['marker_color'], edgecolor=B['marker_edge'], linewidth=0.5, zorder=3)
    max_val = max(merged['hmm_count'].max(), merged['karyoscope_count'].max()) + 5
    ax.plot([0, max_val], [0, max_val], '--', color='#999999', linewidth=1)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ts = B['text_size']
    ax.set_xlabel('HMM repeat count', fontsize=ts, fontweight='bold', labelpad=2.1)
    ax.set_ylabel('KaryoScope repeat count', fontsize=ts, fontweight='bold', labelpad=2.1)
    ax.tick_params(labelsize=ts, length=2.5, width=0.6, pad=1.5)
    ticks = np.arange(0, max_val + 1, 20)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    sx, sy, sw, sh = B['stat_box']
    fx, fy = _to_axes_frac(B['box'], B['axes'], sx + 1.3, sy + 0.9)
    ax.text(fx, fy, f'r = {r:.3f}\nn = {len(merged)}', transform=ax.transAxes, fontsize=ts,
            va='top', ha='left', linespacing=0.95,
            bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='black', linewidth=0.6))
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    path = ctx.panels / 'figS9_B.svg'
    fig.savefig(path, format='svg', transparent=True, metadata={'Date': None})
    plt.close(fig)
    bx, by, bw, bh = B['box']
    ctx.items.append(fc.nest_svg(path, bx, by, bw, bh, prefix='S9B_'))


def build_figS9(ctx):
    ctx.items.append(figure_schematic.svg_fragment(ctx.L['page_origin'], ctx.sans,
                                                   ctx.L.get('schematic_text_size')))
    hmm_scatter(ctx)
    ctx.all_tracks()
    ctx.track_legend()
    ctx.letters()


def build_figS10(ctx):
    ctx.all_tracks()
    ctx.track_legend()
    ctx.letters()


BUILDERS = {'fig4': build_fig4, 'figS9': build_figS9, 'figS10': build_figS10}


def check_haplotypes(name, cfg, layout):
    """Each panel lists at least one haplotype and no more than the layout has track slots for it
    (a panel with fewer leaves its last slots empty)."""
    slots = {}
    for t in layout[name]['tracks']:
        slots[t['panel']] = slots.get(t['panel'], 0) + 1
    for panel, n in slots.items():
        got = len(cfg[name].get(panel, {}).get('haplotypes', []))
        if not 1 <= got <= n:
            raise ValueError(f'{name} panel {panel}: layout has {n} track slots, config lists {got} haplotypes')


def build(name, cfg, layout, out):
    check_haplotypes(name, cfg, layout)
    ctx = FigureContext(name, cfg, layout, out)
    BUILDERS[name](ctx)
    return ctx.write()
