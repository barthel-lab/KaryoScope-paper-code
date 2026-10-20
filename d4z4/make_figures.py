#!/usr/bin/env python3
"""
Regenerate the D4Z4 manuscript figures from the analysis outputs.

    python make_figures.py --config figures_config.yaml --out results/figures/

writes, for each of Figure 4, Supplementary Figure S9 and Supplementary Figure S10:
    <name>.svg   editable vector page (text kept as text)
    <name>.pdf
    <name>.png   at the published PNG's pixel size
plus panels/ (every panel as its own SVG) and figure_data.tsv (every data value drawn).

What is drawn (haplotypes, filters) is set in figures_config.yaml; where it is drawn
(panel positions, font sizes) in figures_layout.yaml, which was read off the published
Illustrator files.

    python make_figures.py --export HG00621#1#CM087849.1 HG00621#2#CM087863.1 --out results/export/HG00621/

instead exports track views of the given haplotypes (sample#hap#accession) in the Supplementary
Figure S10 style, one file each plus side-by-side and stacked versions, and colour legends packed
into several shapes, for use as panels in other figures (see figure_export.py).
"""

import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path

import yaml

BASE_DIR = Path(__file__).resolve().parent


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_paths(cfg, config_path):
    base = Path(config_path).resolve().parent
    for k, v in cfg['inputs'].items():
        if v is not None and not os.path.isabs(v):
            cfg['inputs'][k] = str(base / v)
    if cfg['style'].get('font_dir') and not os.path.isabs(cfg['style']['font_dir']):
        cfg['style']['font_dir'] = str(base / cfg['style']['font_dir'])
    return cfg


def setup_fonts(cfg, workdir):
    """Make the configured font (or a metric-compatible fallback) visible to matplotlib and cairo."""
    style = cfg['style']
    font_dir = style.get('font_dir')
    if font_dir:
        # cairosvg resolves fonts through fontconfig: add font_dir to the default configuration
        conf = Path(workdir) / 'fonts.conf'
        default = os.environ.get('FONTCONFIG_FILE') or ''
        include = f'<include ignore_missing="yes">{default}</include>' if default else \
            '<include ignore_missing="yes">/etc/fonts/fonts.conf</include>'
        prefix_conf = Path(sys.prefix) / 'etc' / 'fonts' / 'fonts.conf'
        if prefix_conf.exists() and not default:
            include = f'<include ignore_missing="yes">{prefix_conf}</include>'
        conf.write_text('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd">'
                        f'<fontconfig>{include}<dir>{font_dir}</dir></fontconfig>\n')
        os.environ['FONTCONFIG_FILE'] = str(conf)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.font_manager as fm
    import matplotlib.pyplot as plt
    if font_dir:
        for f in sorted(Path(font_dir).glob('*.[ot]tf')):
            fm.fontManager.addfont(str(f))
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    plt.rcParams['font.family'] = [style['font_family']] + list(style.get('font_fallbacks', []))
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['svg.hashsalt'] = 'd4z4'  # stable element ids, so reruns give identical SVGs
    plt.rcParams['axes.linewidth'] = 0.6
    plt.rcParams['xtick.major.width'] = 0.6
    plt.rcParams['ytick.major.width'] = 0.6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=str(BASE_DIR / 'figures_config.yaml'))
    ap.add_argument('--layout', default=str(BASE_DIR / 'figures_layout.yaml'))
    ap.add_argument('--out', default=str(BASE_DIR / 'results' / 'figures'))
    ap.add_argument('--font-dir', help='directory of .ttf/.otf fonts (overrides style.font_dir); '
                    'use when Arial/Helvetica are not installed, e.g. Liberation Sans')
    ap.add_argument('--figures', nargs='+', default=['fig4', 'figS9', 'figS10'],
                    choices=['fig4', 'figS9', 'figS10', 'sunburst', 'hprc2'],
                    help='sunburst: the HPRC2 chr4 sunburst, all haplotypes and assembly-QC pass only; '
                         'hprc2: every HPRC2 figure (the sunburst and each entry under hprc2: in the config)')
    ex = ap.add_argument_group('haplotype export (instead of the manuscript figures)')
    ex.add_argument('--export', nargs='+', metavar='CONTIG',
                    help='haplotypes to export as S10-style track views, as sample#hap#accession')
    ex.add_argument('--arrangements', nargs='+', default=['single', 'side_by_side', 'stacked'],
                    choices=['single', 'side_by_side', 'stacked'],
                    help='one file per haplotype, all in one row, and/or one per row (default: all three)')
    ex.add_argument('--legend-columns', nargs='+', default=['wide', '1', '2', '3', '4'],
                    help="legend shapes: 'wide' (one column per group) and/or numbers of columns "
                         '(default: wide 1 2 3 4)')
    ex.add_argument('--legend-features', choices=['drawn', 'all'], default='drawn',
                    help='list only the features drawn in the exported tracks (default), or all')
    ex.add_argument('--scale', choices=['shared', 'own'], default='shared',
                    help='draw all exported haplotypes at one kb scale (default) or each in its own window')
    ex.add_argument('--dpi', type=int, default=600, help='PNG resolution for exports (default: 600)')
    args = ap.parse_args()
    for c in args.legend_columns:
        if c != 'wide' and not (c.isdigit() and int(c) > 0):
            ap.error(f"--legend-columns: {c!r} is neither 'wide' nor a positive number")

    cfg = resolve_paths(load_yaml(args.config), args.config)
    if args.font_dir:
        cfg['style']['font_dir'] = os.path.abspath(args.font_dir)
    layout = load_yaml(args.layout)
    out = Path(args.out)
    (out / 'panels').mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        setup_fonts(cfg, tmp)
        if args.export:
            import figure_export
            tracks, legends = figure_export.export(cfg, layout, args.export, out, args.arrangements,
                                                   args.legend_columns, args.dpi, args.legend_features,
                                                   args.scale)
            for name, paths in tracks.items():
                print(f'{name}: ' + ', '.join(str(p) for p in paths))
            for name, (w, h, white) in legends.items():
                print(f'{name}: {w} x {h} pt, {white:.0%} white space between groups')
            return
        import figure_build
        import figure_data
        summary = figure_data.data_summary(cfg)
        with open(out / 'figure_data.tsv', 'w') as f:
            f.write('key\tvalue\n')
            for k, v in summary.items():
                f.write(f'{k}\t{v}\n')
        for name in args.figures:
            if name == 'sunburst':
                import figure_sunburst
                paths = figure_sunburst.build(cfg, out)
            elif name == 'hprc2':
                import figure_export
                import figure_sunburst
                paths = figure_sunburst.build(cfg, out)
                for fig, spec in cfg.get('hprc2', {}).items():
                    paths += figure_export.composed_figure(fig, spec, cfg, layout, out)
            else:
                paths = figure_build.build(name, cfg, layout, out)
            print(f'{name}: ' + ', '.join(str(p) for p in paths))


if __name__ == '__main__':
    main()
