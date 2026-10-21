# Benchmarks

| directory | manuscript item | source |
|---|---|---|
| `tableS1/` | Supplementary Table 1 (runtime and memory by backend) | cluster job directory `ks_slurm/tableS1_cold/` (in no repository before this one), job 38434006, 2026-10-14 |
| `figS6/` | Extended Data Fig. 6 (thread scaling) | cluster job directory `ks_slurm/threadsweep_cold/`, job 38373939, 2026-10-12; `plot_figS6.py` (the Nature-size version, PR #35) and its two TSVs from `barthel-lab/KaryoScope-paper` `main` `605aa7f`, `figures/fig3/data/` |
| `a1-comparator-runtimes/` | Extended Data Table 1 (comparator tools) | `barthel-lab/KaryoScope-revisions` `b7d2092` |
| `fig2/` | Fig. 2 (concordance with RepeatMasker and CenSat, HG002) | KaryoScope-archive working tree, see below |
| `edfig5/` | Extended Data Fig. 5 (CHM13 self-benchmark) | KaryoScope-archive working tree, see below |

## Table S1 and Extended Data Fig. 6

Both used the same protocol: one job on one Intel Xeon Gold 6342 node (Ice
Lake, 2.80 GHz), points in a fixed random order (`manifest.tsv`; seeds 20261013
and 20261012), and before each point the database, input, binaries and source
evicted from the page cache (`evict.py`, checked with `resident.sh`), then a
cold run and a warm run. The tables report the warm runs; n = 3 per cell.

- `run_tableS1.sh`, `run_sweep.sh`: the SLURM drivers as run. They call
  KaryoScope v2.3.0 (tag commit `124f84a`) with `hks` 0.4.0 (HKS commit
  `04bbc52`) and `get_featureIDs` from v2.3.0. Those two binaries are not
  included; build them from the tagged sources.
- `runs/`: the per-run logs and `runs/index.tsv`. Cluster paths inside the logs
  are replaced by the variables in the repository README.
- `extract_tableS1.py` reads `runs/` and writes `tableS1_runs.tsv` and
  `tableS1_summary.tsv`. Re-running it on the shipped logs reproduces both
  files byte for byte.
- `extract_phases.sh {warm,cold}` reads `runs/` and prints the per-phase
  timings. Its output equals `figS6_phase_timings.{warm,cold}.tsv`.
  `figS6_phase_timings.tsv` (the file the paper repository plots) is the warm
  table.
- `plot_figS6.py --tsv figS6_phase_timings.tsv --out <file>.svg --summary <file>.tsv`
  draws the figure (`environment/figures.yml`). The summary it writes equals
  `figS6_summary.tsv`, and the SVG has the same 116 shapes and 60 texts as
  `figures/fig3/repeat_speedup_icelake.svg` at KaryoScope-paper `605aa7f`, the
  submitted Extended Data Fig. 6.
- `STATUS.md` and `HANDOFF-*.md` are the job notes written at the time.

## Fig. 2 and Extended Data Fig. 5

These R scripts were run from a working tree of the KaryoScope-archive
Snakemake pipeline at commit `7c976e04`, where they were never committed:

| file | sha256 (as copied) | note |
|---|---|---|
| `fig2/figure2_benchmarking.R` | `19f39bd6…34d28346` | rule `figure2_benchmarking` |
| `edfig5/supplementary_CHM13_benchmarking.nature.R` | `b8a05656…d1ead3a0a9bf` | draws the submitted ED Fig. 5: `supplementary_CHM13_benchmarking.R` with Nature styling (Arial, 5–7 pt text, 8 pt lowercase panel letters, 183 × 160 mm, legend and tile-label layout, no overall title). Plotted data identical to the unmodified script |
| `edfig5/supplementary_CHM13_benchmarking.R` | `587d359c…6b49b3bc` | the working-tree script, with the `tandem_repeat` category |
| `edfig5/supplementary_CHM13_benchmarking_old.R` | `66069e78…582c6c87` | the earlier version, without `tandem_repeat`, that drew the first-submission figure |
| `fig2/benchmark.smk` | `b41cf884…1a55aa90` | the working-tree rules file (uncommitted changes), for the rule definitions |
| `fig2/contrast_HG002.py` | `70b2b695…52e1560a` | classifies truth vs prediction per region using the feature hierarchy (`*.contrast.tsv`) |
| `fig2/trim_bed_ends.py` | `0c61a7f8…7192509d` | trims the last k − 1 bases of each sequence from a tiled truth BED |

Inputs are the `*.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv` and
`*.censat.vs.presmoothed.vs.smoothed.contrast.tsv` tables and
`KS_human_CHM13_v2.hierarchy.txt` (in
`vendored/KaryoScope-archive-resources/`). Run with `environment/benchmarks-r.yml`:

```bash
Rscript edfig5/supplementary_CHM13_benchmarking.nature.R --rm <CHM13 repeatmasker contrast.tsv> \
    --cs <CHM13 censat contrast.tsv> --hierarchy <KS_human_CHM13_v2.hierarchy.txt> --outdir out/
Rscript fig2/figure2_benchmarking.R --rm_hap1 ... --rm_hap2 ... --cs_hap1 ... --cs_hap2 ... \
    --hierarchy <KS_human_CHM13_v2.hierarchy.txt> --outdir out/
```

Checked 2026-10-20: the `_old.R` render reproduces every number printed in the
first-submission PDF. The `.nature.R` render plots the same data as the
unmodified `.R` (271 cells); against `_old.R` it adds the `tandem_repeat` bar
(panel a) and row and column (panel c), and every shared cell keeps its printed
label.
The re-rendered Fig. 2 PDF has the same text as the R output kept in the paper
repository; the published Fig. 2 is that output edited in Illustrator. The
per-base concordance numbers themselves are reproduced by
`revision/concordance/ks_concord.sh`.
