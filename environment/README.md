# Environments

Each file is `conda env export --no-builds` from an environment that was built
and then used to run the analysis named below. Each was then rebuilt from the
exported file and the check re-run (2026-10-20, Linux x86_64).

| file | analysis | check run in it |
|---|---|---|
| `centromere.yml` | `centromere/` (Fig. 5a, 5c, 5d) | steps 1-7; 48-row outlier table byte-identical; 4,551 haplotypes, 821 minor |
| `figures.yml` | `vendored/` + `hg002/` (Fig. 1b-e), Fig. 3, `benchmarks/figS6/plot_figS6.py` (ED Fig. 6), the ED Fig. 4 panel scripts | renders compared with the original code and the published files |
| `benchmarks-r.yml` | `benchmarks/fig2/`, `benchmarks/edfig5/` (R) | ED Fig. 5 and Fig. 2 re-rendered |
| `d4z4.yml` | `d4z4/` (Fig. 4, ED Fig. 8, Supp. Fig. 2) | `make_figures.py`; SVGs byte-identical to the HPRC-D4Z4 `77e804b` render |

`centromere.yml` pins the versions that produced the published centromere
counts: Python 3.12.13, numpy 2.4.3, pandas 3.0.1, scipy 1.17.1,
scikit-learn 1.8.0, matplotlib 3.10.8. With other scipy or scikit-learn
versions the minor-haplotype count changes (833 instead of 821 was observed).

Run with `PYTHONNOUSERSITE=1` so that packages in `~/.local` do not shadow the
environment, and for R with `R_LIBS_USER` unset or pointing nowhere.

`figures.yml` installs `karyoplot` from
`git+https://github.com/barthel-lab/KaryoScope-plotlib.git@f59c5fd`.

The revision analyses (`revision/`) and the benchmark drivers call the
KaryoScope command-line tool (v2.3.0) and HKS 0.4.0; install those from their
own repositories (`barthel-lab/KaryoScope`, `jnalanko/HKS` at `04bbc52`).
