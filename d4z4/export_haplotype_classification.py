#!/usr/bin/env python3
"""
Export one row per haplotype with its D4Z4 configuration category and assembly-QC status.

Categories match the waffle charts (create_d4z4_figure.categorize_haplotype); the manuscript-style
category matches the k-sweep counts in the Methods (A-type / B-type / multi-array / degraded / no-array).

Assembly QC (optional) uses results/d4z4_assembly_qc/flag_intervals_by_zone.tsv: a haplotype fails
when an HMM-Flagger Err or NucFlag COLLAPSE/COLLAPSE_VAR interval overlaps its D4Z4 array or terminal
element (Err confined to the last 30 kb of a no-array contig is ignored).

Usage:
    python export_haplotype_classification.py --output results/d4z4_haplotype_classification.tsv
"""
import argparse
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

from create_d4z4_figure import build_plam_lookup, categorize_haplotype

BASE_DIR = Path(__file__).parent
SEVERE = ['Err', 'COLLAPSE', 'COLLAPSE_VAR']
CORE_ZONES = ['array', 'terminal', 'no_array']
TELOMERE_END = 470_000  # window coordinate beyond which a no-array Err is treated as a telomere-end artefact


def manuscript_category(row):
    if row['num_arrays'] == 0:
        return 'no-array'
    if row['terminal_feature'] == 'degraded (no terminal)':
        return 'degraded'
    if row['num_arrays'] > 1:
        return 'multi-array'
    return 'A-type' if row['terminal_feature'] == 'β-satellite terminal' else 'B-type'


def qc_columns(qc_file):
    iv = pd.read_csv(qc_file, sep='\t')
    iv = iv[iv.zone.isin(CORE_ZONES) & ~((iv.zone == 'no_array') & (iv.start >= TELOMERE_END)) & (iv.flag != 'HET')]
    iv['label'] = iv.tool.str.replace('flagger_', '') + ':' + iv.flag  # e.g. hifi:Err, ont:Col, nucflag:MISJOIN
    flags = iv.groupby('haplotype').label.apply(lambda s: ';'.join(sorted(set(s))))
    fail = set(iv[iv.flag.isin(SEVERE)].haplotype)
    return flags, fail


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default=str(BASE_DIR / 'results' / 'd4z4_haplotype_classification.tsv'))
    parser.add_argument('--qc', default=str(BASE_DIR / 'results' / 'd4z4_assembly_qc' / 'flag_intervals_by_zone.tsv'))
    args = parser.parse_args()

    log_file = BASE_DIR / 'results' / f"export_haplotype_classification_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(log_file), logging.StreamHandler()])

    df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    plam_lookup = build_plam_lookup()
    df['category'] = df.apply(lambda r: categorize_haplotype(r, plam_lookup), axis=1)
    df['manuscript_category'] = df.apply(manuscript_category, axis=1)
    df['contig'] = df.haplotype.str.split('_', n=2).str[2]

    plam = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_coordinates.plam.tsv', sep='\t')
    func = set(plam.loc[plam.plam_functional == 'func_Plam', 'contig'])
    nonfunc = set(plam.loc[plam.plam_functional == 'nonfunc_Plam', 'contig'])
    acc = df.contig.str.split('#').str[2]
    df['plam'] = ['functional' if a in func else 'non-functional' if a in nonfunc else 'none' for a in acc]

    cols = ['chromosome', 'sample', 'hap_num', 'contig', 'category', 'manuscript_category', 'num_arrays',
            'num_canonical', 'num_degraded', 'repeat_count_main', 'terminal_feature', 'plam', 'is_truncated']
    if Path(args.qc).exists():
        flags, fail = qc_columns(args.qc)
        df['qc_pass'] = ~df.haplotype.isin(fail)
        df['qc_flags_array_terminal'] = df.haplotype.map(flags).fillna('')
        cols += ['qc_pass', 'qc_flags_array_terminal']
    else:
        logging.warning('QC file not found (%s); writing without QC columns', args.qc)

    out = df[cols].rename(columns={'hap_num': 'hap', 'repeat_count_main': 'repeat_units'})
    out = out.sort_values(['chromosome', 'sample', 'hap'], key=lambda s: s.str.replace('chr', '').str.zfill(3) if s.name == 'chromosome' else s)
    out.to_csv(args.output, sep='\t', index=False)

    logging.info('Wrote %d haplotypes to %s', len(out), args.output)
    logging.info('Categories:\n%s', out.groupby(['chromosome', 'category']).size().to_string())
    if 'qc_pass' in out:
        logging.info('QC pass: %d / %d', out.qc_pass.sum(), len(out))


if __name__ == '__main__':
    main()
