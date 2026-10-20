#!/usr/bin/env python3
"""
Data behind the D4Z4 manuscript figures (Figure 4, Supplementary Figures S9 and S10).

Every number drawn in a data panel is computed here, from the analysis outputs named in
figures_config.yaml. The selection rules (which haplotypes, which filters) are read from the
config so that a change of haplotype list or filter is a one-line edit there.

The classification itself is not reimplemented: categorize_haplotype() is imported from
create_d4z4_figure.py, which is what the published panels were drawn with.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from create_d4z4_figure import categorize_haplotype


def qc_passing(cfg):
    """Contigs (sample#hap#accession) that pass the array-focused assembly QC (qc_pass)."""
    qc = pd.read_csv(cfg['inputs']['qc_classification_tsv'], sep='\t')
    return set(qc.loc[qc['qc_pass'].astype(str) == 'True', 'contig'])


def load_analysis(cfg):
    """The analysis table; with filters: assembly_qc: true, only haplotypes that pass the
    array-focused assembly QC, so every figure panel and number is computed from those."""
    df = pd.read_csv(cfg['inputs']['analysis_tsv'], sep='\t')
    if cfg.get('filters', {}).get('assembly_qc', False):
        passing = qc_passing(cfg)
        df = df[df['haplotype'].str.split('_', n=2).str[2].isin(passing)].reset_index(drop=True)
    return df


def plam_lookup(cfg, analysis_df):
    """Haplotypes with a functional pLAM on chr10 (same rule as build_plam_lookup())."""
    plam_file = Path(cfg['inputs']['plam_tsv'])
    if not plam_file.exists():
        raise FileNotFoundError(f'pLAM table not found: {plam_file}')
    plam_df = pd.read_csv(plam_file, sep='\t')
    chr10_func = plam_df[(plam_df['chromosome'] == 'chr10') &
                         (plam_df['plam_functional'] == 'func_Plam')]
    func_haps = set()
    for _, row in chr10_func.iterrows():
        match = analysis_df[(analysis_df['sample'] == row['sample']) &
                            (analysis_df['chromosome'] == 'chr10') &
                            (analysis_df['hap_num'] == row['hap'])]
        if len(match) > 0:
            func_haps.add(match.iloc[0]['haplotype'])
    return func_haps


def with_config(cfg, analysis_df=None):
    """Analysis table with a 'config' column (structural category per haplotype)."""
    df = load_analysis(cfg) if analysis_df is None else analysis_df.copy()
    lookup = plam_lookup(cfg, df)
    df['config'] = df.apply(lambda row: categorize_haplotype(row, lookup), axis=1)
    return df


def waffle_counts(cfg, df):
    """{chrom: [(category, count), ...]} in display order, plus totals."""
    out = {}
    for chrom, cats in cfg['fig4']['B']['categories'].items():
        cdf = df[df['chromosome'] == chrom]
        out[chrom] = {'n': len(cdf),
                      'counts': [(c, int((cdf['config'] == c).sum())) for c in cats]}
        missing = len(cdf) - sum(n for _, n in out[chrom]['counts'])
        if missing:
            raise ValueError(f'{chrom}: {missing} haplotypes fall outside the configured categories')
    return out


def diplotype_tables(cfg, df):
    """Lower-triangle pairing counts per chromosome (Fig. 4G).

    hap1/hap2 are arbitrary, so each individual contributes one unordered pair.
    Individuals without exactly one haplotype of each are skipped, as in generate_panel_g().
    """
    out = {}
    for chrom, spec in cfg['fig4']['G']['tables'].items():
        cdf = df[df['chromosome'] == chrom]
        exclude = set(spec.get('exclude_categories', []))
        cats = [c for c in spec['categories'] if c in cdf['config'].values and c not in exclude]
        pairs, total = {}, 0
        for sample in cdf['sample'].unique():
            sdf = cdf[cdf['sample'] == sample]
            h1, h2 = sdf[sdf['hap_num'] == 'hap1'], sdf[sdf['hap_num'] == 'hap2']
            if len(h1) != 1 or len(h2) != 1:
                continue
            a, b = sorted([h1.iloc[0]['config'], h2.iloc[0]['config']],
                          key=lambda x: cats.index(x) if x in cats else 999)
            pairs[(a, b)] = pairs.get((a, b), 0) + 1
            total += 1
        n = len(cats)
        matrix = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(i + 1):
                matrix[i, j] = pairs.get((cats[j], cats[i]), 0)
        # Individuals whose pair involves an excluded category are counted in n but not drawn
        out[chrom] = {'categories': cats, 'matrix': matrix, 'n': total}
    return out


def asymmetry_pairs(cfg, df):
    """Shorter vs longer haplotype repeat count per individual (Fig. 4H)."""
    spec = cfg['fig4']['H']
    rows = []
    for chrom in spec['chromosomes']:
        cdf = df[df['chromosome'] == chrom]
        for sample in cdf['sample'].unique():
            sdf = cdf[cdf['sample'] == sample]
            h1, h2 = sdf[sdf['hap_num'] == 'hap1'], sdf[sdf['hap_num'] == 'hap2']
            if len(h1) != 1 or len(h2) != 1:
                continue
            if spec['single_array_only'] and (h1.iloc[0]['num_arrays'] != 1 or
                                              h2.iloc[0]['num_arrays'] != 1):
                continue
            r1, r2 = h1.iloc[0]['repeat_count_main'], h2.iloc[0]['repeat_count_main']
            if pd.isna(r1) or pd.isna(r2) or r1 == 0 or r2 == 0:
                continue
            shorter, longer = sorted([r1, r2])
            rows.append({'sample': sample, 'chromosome': chrom,
                         'shorter': shorter, 'longer': longer, 'fold': longer / shorter})
    pairs = pd.DataFrame(rows)
    # Axis limit: largest repeat_count_main over ALL haplotypes, plus 5 (as published), or over
    # the plotted pairs only (H: axis_from: pairs)
    if spec.get('axis_from', 'all_haplotypes') == 'pairs' and len(pairs):
        axis_max = pairs['longer'].max() + 5
    else:
        axis_max = df['repeat_count_main'].dropna().max() + 5
    return pairs, axis_max


def hmm_comparison(cfg, df):
    """KaryoScope vs HMM repeat-unit counts (Fig. S9B), as in create_d4z4_supp1.generate_panel_b()."""
    from scipy import stats
    hmm_df = pd.read_excel(cfg['inputs']['hmm_counts_xlsx'])
    hmm = {str(r['Contig']): int(r['Count']) for _, r in hmm_df.iterrows()}
    spec = cfg['figS9']['B']
    sel = df[(df['chromosome'] == spec['chromosome']) & (df['num_arrays'] == 1)]
    if spec.get('exclude_truncated', True):
        sel = sel[~sel['is_truncated']]
    merged = []
    for _, row in sel.iterrows():
        contig = row['haplotype'].split('_', 2)[-1]
        if contig in hmm:
            merged.append({'haplotype': row['haplotype'],
                           'karyoscope_count': row['repeat_count_main'],
                           'hmm_count': hmm[contig]})
    merged = pd.DataFrame(merged)
    r, p = stats.pearsonr(merged['karyoscope_count'], merged['hmm_count'])
    return merged, r


def haplotype_id(df, contig):
    """'HG01358#1#CM089053.1' -> 'chr4_hap1_HG01358#1#CM089053.1' (the analysis table's key)."""
    hits = df[df['haplotype'].str.endswith('_' + contig)]
    if len(hits) != 1:
        raise KeyError(f'{contig}: {len(hits)} matches in the analysis table (not analysed, or '
                       'excluded by the assembly-QC filter, filters: assembly_qc); choose another haplotype')
    return hits.iloc[0]['haplotype']


def data_summary(cfg):
    """Every data value drawn, as a flat dict (written next to the figures, used by the check)."""
    df = with_config(cfg)
    s = {}
    for chrom, w in waffle_counts(cfg, df).items():
        s[f'4B.{chrom}.n'] = w['n']
        for c, n in w['counts']:
            s[f'4B.{chrom}.{c}'] = n
    for chrom, t in diplotype_tables(cfg, df).items():
        s[f'4G.{chrom}.n'] = t['n']
        cats, m = t['categories'], t['matrix']
        for i in range(len(cats)):
            for j in range(i + 1):
                s[f'4G.{chrom}.{cats[i]} x {cats[j]}'] = int(m[i, j])
    pairs, axis_max = asymmetry_pairs(cfg, df)
    s['4H.n'] = len(pairs)
    for chrom in cfg['fig4']['H']['chromosomes']:
        s[f'4H.{chrom}.n'] = int((pairs['chromosome'] == chrom).sum())
    s['4H.axis_max'] = float(axis_max)
    for fold in cfg['fig4']['H']['fold_lines']:
        s[f'4H.pairs_above_{fold}x'] = int((pairs['fold'] > fold).sum())
    top = pairs.nlargest(cfg['fig4']['H']['label_top_n'], 'fold')
    for k, (_, row) in enumerate(top.iterrows(), 1):
        s[f'4H.outlier{k}'] = f"{row['sample']} ({row['fold']:.1f}x)"
    merged, r = hmm_comparison(cfg, df)
    s['S9B.n'] = len(merged)
    s['S9B.r'] = round(float(r), 3)
    for fig in ('fig4', 'figS9', 'figS10'):
        omitted = cfg[fig].get('omit_panels', [])
        for panel, spec in cfg[fig].items():
            if not isinstance(spec, dict) or panel in omitted:
                continue
            for k, contig in enumerate(spec.get('haplotypes', []), 1):
                row = df[df['haplotype'] == haplotype_id(df, contig)].iloc[0]
                s[f'{fig}.{panel}.hap{k}'] = (f"{contig} {row['config']} arrays={row['num_arrays']} "
                                              f"units={row['repeat_count_main']}")
    return s
