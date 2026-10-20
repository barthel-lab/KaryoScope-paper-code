#!/usr/bin/env python3
"""
Recompute the manuscript's D4Z4 headline numbers under alternative assembly-QC filters.

Requires results/d4z4_assembly_qc/flag_intervals_by_zone.tsv (run d4z4_assembly_qc.py first).

Filters:
  published      all CM contigs, no QC
  window         no Err/COLLAPSE/COLLAPSE_VAR anywhere in the 500 kb window (the centromere rule);
                 needs the KaryoScope aggregate QC table (--aggregate), skipped if absent
  array-focused  no Err/COLLAPSE/COLLAPSE_VAR overlapping the array or terminal element
                 (Err in the last 30 kb of a no-array contig ignored)  <- recommended
  strict         array-focused + no NucFlag MISJOIN in the array or terminal element

Output: results/d4z4_assembly_qc/scenarios.tsv, array_focused_fail.txt

Usage:
    python d4z4_qc_sensitivity.py [--aggregate path/to/aggregate_qc_v4.manual_curation.tsv]
"""
import argparse
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd
from scipy import stats

BASE_DIR = Path(__file__).parent
QC_DIR = BASE_DIR / 'results' / 'd4z4_assembly_qc'
SEVERE = ['Err', 'COLLAPSE', 'COLLAPSE_VAR']
CORE = ['array', 'terminal', 'no_array']
TELOMERE_END = 470_000


def manuscript_category(r):
    if r.num_arrays == 0:
        return 'no-array'
    if r.terminal_feature == 'degraded (no terminal)':
        return 'degraded'
    if r.num_arrays > 1:
        return 'multi-array'
    return 'A-type' if r.terminal_feature == 'β-satellite terminal' else 'B-type'


def summarize(ana, fail, func, hmm):
    d = ana[~ana.haplotype.isin(fail)]
    c4, c10 = d[d.chromosome == 'chr4'], d[d.chromosome == 'chr10']
    o = {'excluded': len(ana) - len(d), 'chr4 haplotypes': len(c4), 'chr10 haplotypes': len(c10)}
    for name, c in [('chr4', c4), ('chr10', c10)]:
        o[f'{name} beta-sat terminal %'] = round(100 * (c.terminal_feature == 'β-satellite terminal').mean(), 1)
        o[f'{name} no beta-sat terminal %'] = round(100 * (c.terminal_feature == 'No β-satellite terminal').mean(), 1)
    for k in ['A-type', 'B-type', 'multi-array', 'degraded', 'no-array']:
        o[f'{k} (n)'] = int((d.ms == k).sum())
    o['multi-array chr4 / chr10'] = f"{(c4.num_arrays > 1).sum()} / {(c10.num_arrays > 1).sum()}"
    o['chr10 functional pLAM (n)'] = func[(func.chromosome == 'chr10') & func.contig.isin(c10.acc)].contig.nunique()
    f4 = func[(func.chromosome == 'chr4') & func.contig.isin(c4.acc)]
    o['FSHD1-like (chr4 functional pLAM, <=10 units)'] = int((f4.groupby('contig').element_count.min() <= 10).sum())
    o['total repeat units'] = int(d.repeat_count_total.sum())
    pairs = []
    for c in [c4, c10]:
        for _, sd in c.groupby('sample'):
            if len(sd) == 2 and (sd.num_arrays == 1).all():
                r = sorted(sd.repeat_count_main)
                if r[0] > 0:
                    pairs.append(r[1] / r[0])
    o['single-array pairs (n)'] = len(pairs)
    o['pairs >3-fold %'] = round(100 * sum(p > 3 for p in pairs) / len(pairs), 1)
    h = c4[(c4.num_arrays == 1) & (~c4.is_truncated)].merge(hmm, left_on='contig', right_on='Contig')
    o['HMM r (n)'] = f"{stats.pearsonr(h.repeat_count_main, h.Count)[0]:.3f} ({len(h)})"
    return o


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--aggregate', default=str(Path.home() / 'Documents/software/KaryoScope/resources/aggregate_qc_v4.manual_curation.tsv'),
                        help='KaryoScope aggregate QC table (for the window filter)')
    args = parser.parse_args()

    log_file = BASE_DIR / 'results' / f"d4z4_qc_sensitivity_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(log_file), logging.StreamHandler()])

    iv = pd.read_csv(QC_DIR / 'flag_intervals_by_zone.tsv', sep='\t')
    ana = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    ana['ms'] = ana.apply(manuscript_category, axis=1)
    ana['contig'] = ana.haplotype.str.split('_', n=2).str[2]
    ana['acc'] = ana.contig.str.split('#').str[2]
    plam = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_coordinates.plam.tsv', sep='\t')
    func = plam[plam.plam_functional == 'func_Plam']
    hmm = pd.read_excel(BASE_DIR / 'Hailey-notes' / 'HPRC_chr4_counts.xlsx')

    core = iv.zone.isin(CORE) & ~((iv.zone == 'no_array') & (iv.start >= TELOMERE_END))
    array_fail = set(iv[core & iv.flag.isin(SEVERE)].haplotype)
    scenarios = {'published': set(), 'array-focused': array_fail,
                 'strict': array_fail | set(iv[core & (iv.flag == 'MISJOIN')].haplotype)}
    if Path(args.aggregate).exists():
        v4 = pd.read_csv(args.aggregate, sep='\t')
        v4 = v4[(v4.region == 'subtelo_q') & ((v4._n_Err > 0) | (v4._n_COLLAPSE > 0) | (v4._n_COLLAPSE_VAR > 0))]
        scenarios = {'published': set(), 'window': set(ana[ana.contig.isin(v4.contig)].haplotype),
                     **{k: v for k, v in scenarios.items() if k != 'published'}}
    else:
        logging.warning('Aggregate QC table not found (%s); skipping window filter', args.aggregate)

    res = pd.DataFrame({k: summarize(ana, f, func, hmm) for k, f in scenarios.items()})
    res.to_csv(QC_DIR / 'scenarios.tsv', sep='\t')
    pd.Series(sorted(array_fail)).to_csv(QC_DIR / 'array_focused_fail.txt', index=False, header=False)
    logging.info('Scenario summary:\n%s', res.to_string())


if __name__ == '__main__':
    main()
