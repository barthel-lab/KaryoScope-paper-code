#!/usr/bin/env python3
"""Table S1 columns from the runs/ logs, using the same phase definitions as the
published Table S1 (job 35247555):

  assignment  KMC `ran get_featureIDs in Xs`; HKS sum of `hks lookup for '<set>' ... in Xs`
  smoothing   KMC `smoothing pass complete in Xs`; HKS sum of `hks smooth for '<set>' ... in Xs`
  compression sum of `bgzipped ... in Xs`
  end-to-end  /usr/bin/time -v "Elapsed (wall clock)"
  peak RSS    /usr/bin/time -v "Maximum resident set size" (largest single process), GiB

Writes tableS1_runs.tsv (one row per run) and tableS1_summary.tsv
(mean and sd, minutes, by mode x input x backend).
"""
import re
from pathlib import Path

import pandas as pd

W = Path(__file__).resolve().parent


def secs(clock):
    p = [float(x) for x in clock.split(':')]
    return sum(v * 60 ** i for i, v in enumerate(reversed(p)))


def parse(log):
    t = Path(log).read_text()
    if re.search(r'hks lookup for', t):
        a = sum(map(float, re.findall(r"hks lookup for '\w+' wrote .* in ([\d.]+)s", t)))
        s = sum(map(float, re.findall(r"hks smooth for '\w+' wrote .* in ([\d.]+)s", t)))
        nsets = len(re.findall(r"hks smooth for '\w+' wrote", t))
    else:
        a = float(re.search(r'ran get_featureIDs in ([\d.]+)s', t).group(1))
        s = float(re.search(r'smoothing pass complete in ([\d.]+)s', t).group(1))
        nsets = len(set(re.findall(r'bgzipped \S+?\.(\w+)\.smoothed\.bed\.gz', t)))
    c = sum(map(float, re.findall(r'bgzipped .* in ([\d.]+)s', t)))
    e2e = secs(re.search(r'Elapsed \(wall clock\).*: (\S+)', t).group(1))
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', t).group(1)) / 2**20
    exit_status = int(re.search(r'Exit status: (\d+)', t).group(1))
    return dict(assign_s=a, smooth_s=s, compress_s=c, e2e_s=e2e, rss_gib=rss,
                feature_sets=nsets, exit=exit_status)


def main():
    idx = pd.read_csv(W / 'runs' / 'index.tsv', sep='\t')
    runs = pd.concat([idx, pd.DataFrame([parse(W / l) for l in idx.log])], axis=1)
    runs.drop(columns='log').to_csv(W / 'tableS1_runs.tsv', sep='\t', index=False, float_format='%.3f')
    cols = ['assign_s', 'smooth_s', 'compress_s', 'e2e_s']
    g = runs.groupby(['mode', 'input', 'backend'])
    summ = pd.concat([(g[cols].mean() / 60).add_suffix('_min_mean'),
                      (g[cols].std() / 60).add_suffix('_min_sd'),
                      g.rss_gib.agg(['mean', 'std']).add_prefix('rss_gib_'),
                      g.size().rename('n')], axis=1)
    summ.to_csv(W / 'tableS1_summary.tsv', sep='\t', float_format='%.3f')
    print(summ.round(2).to_string())


if __name__ == '__main__':
    main()
