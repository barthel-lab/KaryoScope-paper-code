#!/usr/bin/env python3
"""
Read-based assembly QC of the D4Z4 haplotypes using HPRC Release 2 HMM-Flagger and NucFlag calls.

For each chr4/chr10 haplotype in results/d4z4_array_analysis.tsv this script:
  1. downloads (and caches) the sample's QC BEDs from the public HPRC S3 bucket:
       - HMM-Flagger v1.2.0 conservative predictions, PacBio HiFi and ONT alignments
       - NucFlag v0.3.3 (HiFi), hap1/hap2 or mat/pat naming
  2. places every flagged interval in the terminal 500 kb analysis window (the KaryoScope window,
     [contig_length - 500 kb, contig_length), 0-based half-open)
  3. assigns it to a zone relative to the called D4Z4 arrays:
       proximal | array (first array start to last array end, incl. inter-array gaps)
       | terminal (10 kb distal to the last array, or the beta-satellite terminal if longer) | distal
       | no_array (whole window, for haplotypes without arrays)

Outputs (results/d4z4_assembly_qc/):
  flag_intervals_by_zone.tsv  one row per flagged interval x zone
  contigs.tsv                 one row per haplotype (category, array span, QC coverage)
  qc_sources_manifest.tsv     every QC BED used: sample, tool, source URL, size, SHA-256
  beds/                       downloaded QC BEDs (cache; delete to re-download)

Usage:
    python d4z4_assembly_qc.py [--cache results/d4z4_assembly_qc/beds] [--threads 16]
"""
import argparse
import hashlib
import logging
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import pandas as pd

from create_d4z4_figure import build_plam_lookup, categorize_haplotype

BASE_DIR = Path(__file__).parent
OUT_DIR = BASE_DIR / 'results' / 'd4z4_assembly_qc'
S3 = 'https://s3-us-west-2.amazonaws.com/human-pangenomics'
FLAGGER_PREFIX = 'submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/{s}/hprc_r2/assembly_qc/hmm_flagger/v1.2.0_{tech}/'
NUCFLAG_PREFIX = 'submissions/5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG/{s}/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/'

WINDOW = 500_000
TERMINAL_ZONE = 10_000


def s3_keys(prefix):
    xml = subprocess.run(['curl', '-s', f'{S3}?list-type=2&prefix={prefix}'], capture_output=True, text=True).stdout
    return re.findall(r'<Key>([^<]*)</Key>', xml)


def download(url, dest):
    """Download url to dest (skipped if cached); record the source next to the file for the manifest."""
    src = dest.with_suffix('.bed.source')
    if not dest.exists():
        if subprocess.run(['curl', '-sf', '-o', str(dest), url]).returncode != 0:
            return False
        src.write_text(url + '\n')
    return True


def fetch(sample, cache):
    """Download the Flagger (HiFi, ONT) and NucFlag BEDs for one sample; return list of local paths.

    Files are resolved by listing the sample's S3 prefix, because file names differ between samples
    (e.g. NucFlag assembly versions v1.0.1 / v1.1.0 / v1, and hap1/hap2 vs mat/pat for trio samples).
    """
    paths = []
    for tech in ['hifi', 'ont']:
        dest = cache / f'{sample}.flagger_{tech}.bed'
        if not dest.exists():
            keys = [k for k in s3_keys(FLAGGER_PREFIX.format(s=sample, tech=tech)) if k.endswith('_prediction.conservative.bed')]
            if not keys or not download(f'{S3}/{keys[0]}', dest):
                logging.warning('%s: no HMM-Flagger %s prediction', sample, tech)
                continue
        paths.append(dest)
    existing = sorted(cache.glob(f'{sample}.nucflag_*.bed'))
    if not existing:
        keys = [k for k in s3_keys(NUCFLAG_PREFIX.format(s=sample)) if k.endswith('.nucflag.bed')]
        # one file per haplotype; the first in sort order is the earliest assembly version (v1.0.1 before v1.0.2)
        by_hap = {}
        for k in sorted(keys):
            m = re.search(rf'{sample}_(hap1|hap2|mat|pat)_', k)
            if m:
                by_hap.setdefault(m.group(1), k)
        for hap, k in by_hap.items():
            dest = cache / f'{sample}.nucflag_{hap}.bed'
            if download(f'{S3}/{k}', dest):
                existing.append(dest)
        if not by_hap:
            logging.warning('%s: no NucFlag v0.3.3 BED', sample)
    return paths + existing


def write_manifest(files, cache, out):
    """One row per QC BED used: sample, tool, local file, source URL, size and SHA-256."""
    rows = []
    for sample, paths in sorted(files.items()):
        for p in paths:
            src = p.with_suffix('.bed.source')
            rows.append(dict(sample=sample, tool=p.name.split('.')[1], file=p.name,
                             url=src.read_text().strip() if src.exists() else '',
                             bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    missing = sum(r['url'] == '' for r in rows)
    if missing:
        logging.warning('%d cached files have no recorded source URL (downloaded before manifests); '
                        'delete %s and rerun to re-download with provenance', missing, cache)


def read_bed(path):
    tool = 'nucflag' if '.nucflag_' in path.name else path.name.split('.')[1]
    d = pd.read_csv(path, sep='\t', header=None, comment='t', usecols=[0, 1, 2, 3], names=['contig', 'start', 'end', 'flag'])
    d['tool'] = tool
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cache', default=str(OUT_DIR / 'beds'), help='Directory for downloaded QC BEDs')
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    log_file = BASE_DIR / 'results' / f"d4z4_assembly_qc_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(log_file), logging.StreamHandler()])

    ana = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    plam_lookup = build_plam_lookup()
    ana['cat'] = ana.apply(lambda r: categorize_haplotype(r, plam_lookup), axis=1)
    ana['contig'] = ana.haplotype.str.split('_', n=2).str[2]
    coords = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_coordinates.tsv', sep='\t')

    samples = sorted(ana['sample'].unique())
    logging.info('Fetching QC BEDs for %d samples into %s', len(samples), cache)
    with ThreadPoolExecutor(args.threads) as ex:
        files = dict(zip(samples, ex.map(lambda s: fetch(s, cache), samples)))

    write_manifest(files, cache, OUT_DIR / 'qc_sources_manifest.tsv')

    rows, contig_rows = [], []
    for sample, sdf in ana.groupby('sample'):
        beds = pd.concat([read_bed(p) for p in files[sample]]) if files[sample] else pd.DataFrame(columns=['contig', 'start', 'end', 'flag', 'tool'])
        for _, h in sdf.iterrows():
            b = beds[beds.contig == h.contig]
            flagger = b[b.tool.str.startswith('flagger')]
            length = flagger.end.max() if len(flagger) else None  # Flagger tiles the whole contig
            acc = h.contig.split('#')[2]
            a = coords[(coords.chromosome == h.chromosome) & (coords.contig == acc)]
            if len(a):
                a0, a1 = a.array_start.min(), a.array_end.max()
                t1 = max(a1 + TERMINAL_ZONE, a.bsat_end.max() if a.bsat_end.notna().any() else 0)
                zones = [('proximal', 0, a0), ('array', a0, a1), ('terminal', a1, t1), ('distal', t1, WINDOW)]
            else:
                a0 = a1 = None
                zones = [('no_array', 0, WINDOW)]
            contig_rows.append(dict(
                haplotype=h.haplotype, contig=h.contig, chromosome=h.chromosome, cat=h['cat'],
                num_arrays=h.num_arrays, repeat_units=h.repeat_count_main, terminal=h.terminal_feature,
                contig_len=length, array_start=a0, array_end=a1,
                has_flagger_hifi=(b.tool == 'flagger_hifi').any(), has_flagger_ont=(b.tool == 'flagger_ont').any(),
                has_nucflag=(b.tool == 'nucflag').any()))
            if length is None:
                logging.warning('%s: no Flagger intervals, contig length unknown; skipped', h.haplotype)
                continue
            w0 = length - WINDOW
            b = b[(b.end > w0) & (b.flag != 'Hap')]
            for _, x in b.iterrows():
                s, e = max(x.start, w0) - w0, x.end - w0
                for zone, z0, z1 in zones:
                    overlap = min(e, z1) - max(s, z0)
                    if overlap > 0:
                        rows.append(dict(haplotype=h.haplotype, chromosome=h.chromosome, cat=h['cat'], tool=x.tool,
                                         flag=x.flag, start=s, end=e, zone=zone, overlap_bp=overlap, dist_to_end=WINDOW - e))

    iv = pd.DataFrame(rows)
    cg = pd.DataFrame(contig_rows)
    iv.to_csv(OUT_DIR / 'flag_intervals_by_zone.tsv', sep='\t', index=False)
    cg.to_csv(OUT_DIR / 'contigs.tsv', sep='\t', index=False)
    logging.info('Haplotypes: %d (Flagger HiFi %d, Flagger ONT %d, NucFlag %d)', len(cg),
                 cg.has_flagger_hifi.sum(), cg.has_flagger_ont.sum(), cg.has_nucflag.sum())
    logging.info('Haplotypes with each flag, by zone:\n%s',
                 iv.groupby(['tool', 'flag', 'zone']).haplotype.nunique().unstack(fill_value=0).to_string())


if __name__ == '__main__':
    main()
