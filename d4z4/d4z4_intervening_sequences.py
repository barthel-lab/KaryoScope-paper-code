#!/usr/bin/env python3
"""
Intervening sequences between the D4Z4 arrays of multi-array and mixed haplotypes.

An intervening sequence is the stretch of assembly between two adjacent arrays of one haplotype
(from the end of one array to the start of the next, in results/d4z4_array_coordinates.tsv).
For every one this script:
  1. extracts the sequence from the assembly FASTA and counts its N bases (an N run means the
     measured length is an assembly gap, not sequence)
  2. records its KaryoScope label composition (region, gene, repeat, subtelomeric feature sets)
  3. aligns it with minimap2 (asm20) to every other intervening sequence, and to 20 kb of D4Z4
     array sequence from one long single array per chromosome (5 kb in from the array start)
  4. groups the sequences of haplotypes that pass assembly QC and contain no N into length
     classes (sorted lengths per chromosome; a new class starts where consecutive lengths differ
     by more than --class-gap bp) and summarises each class

Inputs:
  results/d4z4_array_coordinates.tsv      arrays per haplotype (extract_array_coordinates.py)
  results/d4z4_haplotype_classification.tsv  category and qc_pass (export_haplotype_classification.py)
  data/pangenome.{chr}.subtelo_q.KS_human_CHM13.presmoothed.{set}.bed  KaryoScope annotation
  --assembly-dir  HPRC Release 2 assemblies, bgzipped with .fai and .gzi
                  (<sample>_<hap>_hprc_r2_*.fa.gz; contigs named sample#hap#accession)

All coordinates are window-relative, 0-based half-open, in the terminal 500 kb of the contig.

Outputs (results/d4z4_intervening_sequences/):
  intervening_sequences.tsv  one row per intervening sequence
  intervening_pairs.tsv      best alignment per pair of intervening sequences (each pair once)
  length_classes.tsv         one row per length class (QC-passing, gapless sequences only; the
                             example_ks_* columns are the composition of its shortest member)
  intervening_sequences.fa, d4z4_reference.fa

Requires samtools and minimap2 on PATH (or --samtools / --minimap2).

Usage:
    python d4z4_intervening_sequences.py [--assembly-dir DIR] [--class-gap 1000]
"""
import os as _ks_os
KS_REFERENCES = _ks_os.environ.get("KS_REFERENCES", "$KS_REFERENCES")  # see the repository README, Paths
import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

BASE_DIR = Path(__file__).parent
RESULTS = BASE_DIR / 'results'
OUT_DIR = RESULTS / 'd4z4_intervening_sequences'
WINDOW = 500_000
FEATURE_SETS = ['region', 'gene', 'repeat', 'subtelomeric']
PAIR_MIN_ALIGNED = 0.9    # a pair counts as aligned when the best hit covers >90% of both sequences


def contig_index(assembly_dir):
    """contig name -> (FASTA path, contig length), from every .fa.gz.fai in the directory."""
    idx = {}
    for fai in Path(assembly_dir).glob('*.fa.gz.fai'):
        for line in fai.read_text().splitlines():
            name, length = line.split('\t')[:2]
            idx[name] = (fai.with_suffix(''), int(length))
    return idx


def fetch(samtools, idx, contig, start, end):
    """Window coordinates [start, end) of contig -> sequence."""
    fasta, length = idx[contig]
    w0 = length - WINDOW
    out = subprocess.run([samtools, 'faidx', str(fasta), f'{contig}:{w0 + start + 1}-{w0 + end}'],
                         capture_output=True, text=True, check=True).stdout
    return ''.join(out.splitlines()[1:])


def composition(beds, chrom, hap, start, end):
    """Per feature set: labels covering >=5% of [start, end), as 'label:NN%;...'."""
    out = {}
    for s in FEATURE_SETS:
        b = beds[(chrom, s)]
        b = b[(b.c == hap) & (b.e > start) & (b.s < end)]
        bp = defaultdict(int)
        for x in b.itertuples():
            bp[x.l.replace('_specific', '')] += min(x.e, end) - max(x.s, start)
        tot = end - start
        out[f'ks_{s}'] = ';'.join(f'{k}:{100 * v / tot:.0f}%' for k, v in sorted(bp.items(), key=lambda kv: -kv[1])
                                  if v / tot >= 0.05)
    return out


def paf(minimap2, args):
    out = subprocess.run([minimap2, *args], capture_output=True, text=True, check=True).stdout
    cols = ['q', 'ql', 'qs', 'qe', 'st', 't', 'tl', 'ts', 'te', 'm', 'al', 'mq']
    df = pd.DataFrame([line.split('\t')[:12] for line in out.splitlines()], columns=cols)
    return df.astype({c: int for c in ['ql', 'qs', 'qe', 'tl', 'ts', 'te', 'm', 'al']})


def covered_fraction(g):
    """Fraction of the query covered by the union of its alignments."""
    tot, cur = 0, None
    for s, e in sorted(zip(g.qs, g.qe)):
        if cur is None or s > cur[1]:
            if cur:
                tot += cur[1] - cur[0]
            cur = [s, e]
        else:
            cur[1] = max(cur[1], e)
    tot += cur[1] - cur[0]
    return tot / g.ql.iloc[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--assembly-dir', default=(KS_REFERENCES + '/pangenome/release2'))
    ap.add_argument('--results', default=str(RESULTS))
    ap.add_argument('--data', default=str(BASE_DIR / 'data'))
    ap.add_argument('--out', default=str(OUT_DIR))
    ap.add_argument('--class-gap', type=int, default=1000)
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--minimap2', default='minimap2')
    args = ap.parse_args()
    results, out = Path(args.results), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    co = pd.read_csv(results / 'd4z4_array_coordinates.tsv', sep='\t')
    cls = pd.read_csv(results / 'd4z4_haplotype_classification.tsv', sep='\t')
    cls['acc'] = cls.contig.str.split('#').str[2]
    cls['qc_pass'] = cls.qc_pass.astype(str) == 'True'
    beds = {(c, s): pd.read_csv(Path(args.data) / f'pangenome.{c}.subtelo_q.KS_human_CHM13.presmoothed.{s}.bed',
                                sep='\t', header=None, names=['c', 's', 'e', 'l'])
            for c in ('chr4', 'chr10') for s in FEATURE_SETS}
    idx = contig_index(args.assembly_dir)

    # 1-2. every intervening sequence of a multi-array or mixed haplotype
    rows = []
    with open(out / 'intervening_sequences.fa', 'w') as fa:
        for h in cls[cls.category.isin(['Multi-array', 'Mixed'])].itertuples():
            a = co[(co.chromosome == h.chromosome) & (co.contig == h.acc)].sort_values('array_start')
            hap = f'{h.chromosome}_{h.hap}_{h.contig}'
            for i in range(len(a) - 1):
                left, right = a.iloc[i], a.iloc[i + 1]
                start, end = int(left.array_end), int(right.array_start)
                sid = f'{h.chromosome}|{h.contig}|{i + 1}'
                seq = fetch(args.samtools, idx, h.contig, start, end)
                fa.write(f'>{sid}\n{seq}\n')
                rows.append(dict(id=sid, chromosome=h.chromosome, sample=h.sample, hap=h.hap, contig=h.contig,
                                 category=h.category, qc_pass=h.qc_pass, index=i + 1, start=start, end=end,
                                 length=end - start, N_bp=seq.upper().count('N'),
                                 left_array=f'{left.array_type}/{left.element_count}u',
                                 right_array=f'{right.array_type}/{right.element_count}u',
                                 **composition(beds, h.chromosome, hap, start, end)))
    iv = pd.DataFrame(rows)

    # 3a. D4Z4 reference: 20 kb from a long single A-type array that passes QC, one per chromosome
    with open(out / 'd4z4_reference.fa', 'w') as fa:
        for chrom in ('chr4', 'chr10'):
            h = cls[(cls.chromosome == chrom) & (cls.category == 'A-type (C+)') & cls.qc_pass
                    & (cls.repeat_units >= 25)].iloc[0]
            a = co[(co.chromosome == chrom) & (co.contig == h.acc)].iloc[0]
            fa.write(f'>{chrom}_D4Z4_{h.contig}\n'
                     f'{fetch(args.samtools, idx, h.contig, int(a.array_start) + 5000, int(a.array_start) + 25000)}\n')
    ref = paf(args.minimap2, ['-c', '-x', 'asm20', '-N', '50', '-p', '0.1',
                              str(out / 'd4z4_reference.fa'), str(out / 'intervening_sequences.fa')])
    frac = {q: covered_fraction(g) for q, g in ref.groupby('q')}
    iv['frac_aligned_to_D4Z4'] = iv.id.map(lambda q: round(frac.get(q, 0.0), 3))

    # 3b. all against all (-X reports each pair once): best alignment per pair
    av = paf(args.minimap2, ['-c', '-x', 'asm20', '-X', '-N', '50',
                             str(out / 'intervening_sequences.fa'), str(out / 'intervening_sequences.fa')])
    av = av[av.q != av.t]
    av = av.assign(pair=[tuple(sorted(x)) for x in zip(av.q, av.t)])
    best = av.sort_values('al', ascending=False).drop_duplicates('pair')
    pairs = pd.DataFrame(dict(q=best.q, t=best.t, q_len=best.ql, t_len=best.tl,
                              identity=(best.m / best.al).round(4),
                              q_aligned_frac=((best.qe - best.qs) / best.ql).round(3),
                              t_aligned_frac=((best.te - best.ts) / best.tl).round(3)))

    # 4. length classes over QC-passing, gapless intervening sequences
    ok = iv[iv.qc_pass & (iv.N_bp == 0)].sort_values(['chromosome', 'length'])
    classes = []
    for chrom, g in ok.groupby('chromosome'):
        members, prev = [], None
        for r in g.itertuples():
            if prev is not None and r.length - prev > args.class_gap:
                classes.append((chrom, members))
                members = []
            members.append(r)
            prev = r.length
        classes.append((chrom, members))
    crows = []
    for chrom, members in classes:
        ids = {m.id for m in members}
        p = pairs[pairs.q.isin(ids) & pairs.t.isin(ids)]
        aligned = p[(p.q_aligned_frac > PAIR_MIN_ALIGNED) & (p.t_aligned_frac > PAIR_MIN_ALIGNED)]
        n = len(ids)
        crows.append(dict(chromosome=chrom, n_sequences=n,
                          n_haplotypes=len({m.contig for m in members}),
                          length_min=min(m.length for m in members), length_max=max(m.length for m in members),
                          categories=';'.join(sorted({m.category for m in members})),
                          pairs=n * (n - 1) // 2, pairs_aligned=len(aligned),
                          min_identity=aligned.identity.min() if len(aligned) else None,
                          max_frac_aligned_to_D4Z4=max(m.frac_aligned_to_D4Z4 for m in members),
                          haplotypes=';'.join(sorted({f'{m.sample}#{m.hap[-1]}' for m in members})),
                          example_ks_region=members[0].ks_region, example_ks_repeat=members[0].ks_repeat))

    iv.to_csv(out / 'intervening_sequences.tsv', sep='\t', index=False)
    pairs.to_csv(out / 'intervening_pairs.tsv', sep='\t', index=False)
    lc = pd.DataFrame(crows)
    lc.to_csv(out / 'length_classes.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 70)
    print(lc.drop(columns=['example_ks_region', 'example_ks_repeat']).to_string(index=False))


if __name__ == '__main__':
    main()
