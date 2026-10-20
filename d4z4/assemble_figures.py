#!/usr/bin/env python3
"""
Assemble individual figure panels into composite letter-format pages (8.5x11 portrait).

Uses rsvg-convert for SVG→PNG conversion and matplotlib for composition.

Outputs:
  results/karyoscope_companion/figure5_composite.{png,pdf}
  results/karyoscope_companion/supp_figure1_composite.{png,pdf}
  results/karyoscope_companion/supp_figure2_composite.{png,pdf}
"""

import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results' / 'karyoscope_companion'

plt.rcParams['font.family'] = 'Arial'

PAGE_W, PAGE_H = 8.5, 11
DPI = 300


def svg_to_img(svg_path, width=2400):
    """Convert SVG to numpy image array via rsvg-convert."""
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        subprocess.run(
            ['rsvg-convert', '-w', str(width), '-o', tmp_path, str(svg_path)],
            check=True, capture_output=True
        )
        return mpimg.imread(tmp_path)
    finally:
        os.unlink(tmp_path)


def load_img(path):
    """Load PNG or convert SVG to image array."""
    path = Path(path)
    if path.suffix == '.svg':
        return svg_to_img(path)
    return mpimg.imread(str(path))


def add_label(ax, label, x=-0.02, y=1.02):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='right',
            fontfamily='Arial')


def save_fig(fig, name):
    out = OUTPUT_DIR / name
    fig.savefig(str(out) + '.png', dpi=DPI, facecolor='white', bbox_inches='tight')
    fig.savefig(str(out) + '.pdf', dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}.png/.pdf")


def assemble_figure5():
    fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=DPI)
    gs = GridSpec(5, 2, figure=fig,
                  height_ratios=[1.2, 1.0, 0.8, 0.8, 1.0],
                  hspace=0.15, wspace=0.08,
                  left=0.04, right=0.98, top=0.98, bottom=0.02)

    # (a) Schematic track - full width
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_a_schematic.svg'))
    ax.axis('off'); add_label(ax, 'a')

    # (b) Ridge - left
    ax = fig.add_subplot(gs[1, 0])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_b_ridge.png'))
    ax.axis('off'); add_label(ax, 'b')

    # (b') Waffle chart - right
    ax = fig.add_subplot(gs[1, 1])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_b_waffle.png'))
    ax.axis('off'); add_label(ax, "b'")

    # (c) Multi-array - full width
    ax = fig.add_subplot(gs[2, :])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_c_multi_array.svg'))
    ax.axis('off'); add_label(ax, 'c')

    # (d) Degraded/missing - full width
    ax = fig.add_subplot(gs[3, :])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_d_degraded_missing.svg'))
    ax.axis('off'); add_label(ax, 'd')

    # (e) Asymmetry - full width
    ax = fig.add_subplot(gs[4, :])
    ax.imshow(load_img(OUTPUT_DIR / 'fig5_panel_e_asymmetry.png'))
    ax.axis('off'); add_label(ax, 'e')

    save_fig(fig, 'figure5_composite')


def assemble_supp_figure1():
    fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=DPI)
    gs = GridSpec(3, 2, figure=fig,
                  height_ratios=[1.3, 1.0, 0.7],
                  hspace=0.15, wspace=0.08,
                  left=0.04, right=0.98, top=0.98, bottom=0.02)

    # (a) Representative - full width
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(load_img(OUTPUT_DIR / 'supp1_panel_a_representative.svg'))
    ax.axis('off'); add_label(ax, 'a')

    # (b) HMM correlation - left
    ax = fig.add_subplot(gs[1, 0])
    ax.imshow(load_img(OUTPUT_DIR / 'supp1_panel_b_hmm_correlation.png'))
    ax.axis('off'); add_label(ax, 'b')

    # (c) Pipeline - right
    ax = fig.add_subplot(gs[1, 1])
    png = OUTPUT_DIR / 'supp1_panel_c_pipeline.png'
    svg = OUTPUT_DIR / 'supp1_panel_c_pipeline.svg'
    ax.imshow(load_img(png if png.exists() else svg))
    ax.axis('off'); add_label(ax, 'c')

    # (d) Placeholder
    ax = fig.add_subplot(gs[2, :])
    ax.text(0.5, 0.5, '(d) Classification decision tree [to be added]',
            transform=ax.transAxes, ha='center', va='center',
            fontsize=11, color='#999999', style='italic')
    ax.axis('off'); add_label(ax, 'd')

    save_fig(fig, 'supp_figure1_composite')


def assemble_supp_figure2():
    rows = [
        ('supp2_row1_multi_array_chr4.svg', 'Multi-array (chr4)'),
        ('supp2_row2_multi_array_chr10.svg', 'Multi-array (chr10)'),
        ('supp2_row3_mixed_chr4.svg', 'Mixed (chr4)'),
        ('supp2_row4_degraded_only_chr10.svg', 'Degraded-only (chr10)'),
        ('supp2_row5_functional_plam_chr10.svg', 'Func. pLAM (chr10)'),
    ]
    labels = 'abcde'

    fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=DPI)
    gs = GridSpec(len(rows), 1, figure=fig,
                  hspace=0.06,
                  left=0.04, right=0.95, top=0.98, bottom=0.02)

    for i, (fname, category) in enumerate(rows):
        ax = fig.add_subplot(gs[i, 0])
        svg_path = OUTPUT_DIR / fname
        if svg_path.exists():
            ax.imshow(load_img(svg_path))
        ax.axis('off')
        add_label(ax, labels[i])
        ax.text(1.01, 0.5, category, transform=ax.transAxes,
                fontsize=7, va='center', ha='left', rotation=270,
                fontfamily='Arial', color='#545454')

    save_fig(fig, 'supp_figure2_composite')


def main():
    print("Assembling composite figures (letter format, portrait)")
    print("=" * 55)
    print("\nFigure 5:")
    assemble_figure5()
    print("\nSupplementary Figure 1:")
    assemble_supp_figure1()
    print("\nSupplementary Figure 2:")
    assemble_supp_figure2()
    print(f"\nAll composites in {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
