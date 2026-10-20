#!/usr/bin/env python3
"""
Track-view panels for the D4Z4 manuscript figures.

Draws one haplotype per call with plot_d4z4_arrays.plot_single_haplotype(), exactly as
plot_d4z4_arrays.py does from the command line with --adaptive-window, but without the
legend and timestamp footer, and with the manuscript's fonts and array-label layout.
"""

import os

import drawsvg as draw
import pandas as pd

import plot_d4z4_arrays as pda

ANNOT_TYPES = ('region', 'gene', 'repeat', 'subtelomeric')


def bed_path(bed_dir, chrom, annot_type):
    """The KaryoScope BED for one annotation type (same name as plot_d4z4_arrays.get_haplotype_data())."""
    return os.path.join(bed_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.{annot_type}.bed')


def check_beds(bed_dir, chrom):
    """Fail loudly if a BED is missing or unreadable.

    plot_d4z4_arrays.load_bed_file() returns no rows for a missing file, so a broken symlink in
    data/ (the BEDs are links into the cluster project directory) used to give empty tracks
    without any error.
    """
    for annot_type in ANNOT_TYPES:
        path = bed_path(bed_dir, chrom, annot_type)
        if os.path.islink(path) and not os.path.exists(path):
            raise FileNotFoundError(f'{path} is a broken symlink (target {os.readlink(path)} is not reachable)')
        if not os.path.isfile(path):
            raise FileNotFoundError(f'{path} not found')
        with open(path):   # unreadable (permissions, unmounted share) raises here
            pass


def plam_data_for(plam_df, hap):
    """pLAM coordinates per array for one haplotype (same lookup as plot_d4z4_arrays.main())."""
    chrom, hap_num, rest = hap.split('_', 2)
    sample = rest.split('#')[0]
    rows = plam_df[(plam_df['chromosome'] == chrom) & (plam_df['sample'] == sample) &
                   (plam_df['hap'] == hap_num)]
    if len(rows) == 0:
        return None
    out = {}
    for _, row in rows.iterrows():
        if pd.notna(row.get('plam_start')):
            out[int(row['array_index'])] = {'plam_start': row['plam_start'],
                                            'plam_end': row['plam_end'],
                                            'plam_functional': row['plam_functional']}
    return out


class TrackRenderer:
    def __init__(self, cfg, analysis_df, style):
        self.cfg = cfg
        self.df = analysis_df
        self.plam_df = pd.read_csv(cfg['inputs']['plam_tsv'], sep='\t')
        self.color_maps = pda.load_all_colors(cfg['inputs']['colors_dir'])
        for annot_type in ANNOT_TYPES:
            if not self.color_maps.get(annot_type):
                raise FileNotFoundError(f"no {annot_type} colours in {cfg['inputs']['colors_dir']}")
        self.style = style
        self.used_colors = {}

    def window_size(self, hap):
        """The adaptive window plot_single_haplotype() would draw for this haplotype, in bp."""
        chrom = hap.split('_')[0]
        check_beds(self.cfg['inputs']['bed_dir'], chrom)
        hap_data = pda.get_haplotype_data(hap, self.cfg['inputs']['bed_dir'], chrom)
        raw_max_pos = max(e for t in ANNOT_TYPES for _, e, _ in hap_data[t])
        return pda.adaptive_window_size(raw_max_pos, self.df[self.df['haplotype'] == hap].iloc[0])

    def draw(self, hap, out_svg, width, left_margin, right_margin, height_ratio):
        """Write one haplotype's track view to out_svg; returns the drawing's (width, height)."""
        chrom = hap.split('_')[0]
        check_beds(self.cfg['inputs']['bed_dir'], chrom)
        hap_data = pda.get_haplotype_data(hap, self.cfg['inputs']['bed_dir'], chrom)
        # every analysed haplotype has rows in all four BEDs (582/582 in HPRC Release 2)
        empty = [t for t in ANNOT_TYPES if not hap_data[t]]
        if empty:
            raise ValueError(f'{hap}: no rows in the {", ".join(empty)} BED(s) in '
                             f'{self.cfg["inputs"]["bed_dir"]}; the BED does not match the analysis table')
        row = self.df[self.df['haplotype'] == hap].iloc[0]
        plot_width = width - left_margin - right_margin
        d = draw.Drawing(width, 1000)
        # title baseline: the title font size, as plot_single_haplotype() computes it
        top = self.style.get('title_font') or max(9, int(int(15 * height_ratio) * 0.85))
        end_y = pda.plot_single_haplotype(
            d, hap_data, hap, row, top, left_margin, plot_width, '#000000', 'white',
            self.color_maps, self.used_colors, adaptive_window=True,
            plam_data=plam_data_for(self.plam_df, hap), height_ratio=height_ratio,
            style=self.style)
        height = end_y - 20 + self.style.get('bottom_pad', 0)
        d.height = height
        d.view_box = (0, 0, width, height)
        d.save_svg(out_svg)
        return width, height
