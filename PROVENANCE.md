# Provenance

Where each directory came from, and which software versions produced the
published results. The changes made to each directory are listed in its README.

## Software versions in detail

- **KaryoScope v1.1.0** (tag `4c09ce5`, doi:10.5281/zenodo.20657817): the
  annotations behind Figs. 1–2 (HG002 v1.1, `KS_human_CHM13_v2`) and the D4Z4
  and centromere inputs (HPRC Release 2, `KS_human_CHM13`, presmoothed). These
  were first made with the KaryoScope-archive Snakemake pipeline (`7c976e04`
  for HG002, `0293e85` for HPRC), and v1.1.0 reproduces them.
- **Binning of the D4Z4 and centromere windows** (100 bp and 100 kb):
  `vendored/KaryoScope-archive/bin_features.py` (KaryoScope-archive `4434885`).
  KaryoScope's `bin` command uses a different labelling rule, so these windows
  need `bin_features.py`.
- **KaryoScope v2.3.0** (tag `124f84a`, doi:10.5281/zenodo.21852728) with
  **hks 0.4.0** (HKS `04bbc52`): the benchmarks (Supp. Table 1, ED Fig. 6), the
  tool comparison (ED Table 1), and the *k*-length (A2), divergence (A5), and
  cross-genome (A3) analyses.
  - Supp. Table 1 and ED Fig. 6 ran on 2026-10-12 and 2026-10-13 from a clean
    worktree of the v2.3.0 tag.
  - The other analyses ran between 2026-08-11 and 2026-09-03 from a checkout
    seven or eight commits past the tag (`052055a`, `501d5fd`, which report
    "2.3.0"). Those commits change only the `karyotype`, `scaffold`, and
    haplotype-label code, not `annotate`, `bin`, or `build`.
  - Exception: the A2 databases at *k* = 11–61 were built on 2026-08-03, before
    the v2.3.0 tag, with `karyoscope build` from `9ccba0b` (v2.1.0 + 3 commits)
    and hks `c5fecf8`. The annotations at each *k* used the v2.3.0 code with hks
    `04bbc52`. Those databases were not kept (see
    `revision/a2-k-sensitivity/step3/README.md`).
- **Robertsonian figures** (Fig. 3, ED Fig. 7): the assemblies were first
  annotated, scaffolded, and binned with the KaryoScope-archive pipeline (commit
  not recorded). KaryoScope v1.1.0 (annotate and scaffold) with `karyoscope bin`
  reproduces the submitted Fig. 3 on the tested sample (the Fig. 3a assembly),
  before the 4-unit headroom edit applied by `robs/scripts/fig3_headroom.py`.
  Binning with v1.1.0 or v2.3.0 gives the same Fig. 3 on all three samples (see
  `robs/README.md`).
- **KaryoScope-ISCN analysis modules** used by the A2 and A6 Robertsonian
  analyses: commit `5507369`, vendored in `vendored/KaryoScope-ISCN/src/`
  (`PROVENANCE.md` there).

## Sources

| directory | repository | commit | notes |
|---|---|---|---|
| `d4z4/` | barthel-lab/HPRC-D4Z4 | `77e804b`, figure code `52ed8cc` | heads of PR #1 (`figures/reproducible-layout`) and PR #2 (`nature/ed-spec`) |
| `robs/` | barthel-lab/KaryoScope-ROBs | `1188850` | head of PR #2 (`nature/ed-spec`, on `fig5/assemble-panels`), with the configuration replaced by a template |
| `centromere/scripts/` | barthel-lab/KaryoScope-analysis | `90e2874` | branch `pangenome_structure_V2`, `scripts/` only |
| `centromere/qc/` | KaryoScope-archive `0bf3a7d` (`workflow/rules/pangenome.smk`, `get_regions.py`) and two files never committed (`aggregate_master_qc.py`, `aggregate_qc_v2.awk`) | n/a | originals in `centromere/qc/original/`. The runnable scripts were written for this repository and reproduce the original tables |
| `hg002/` | barthel-lab/KaryoScope-HG002 | `27885ba` | `conf/`, `scripts/check_data.sh`, `requirements.txt` |
| `hg002/scripts/make_centromeres_v1coords.py`, `hg002/data/` | written for this repository | n/a | rebuilds the Fig. 1b centromere window file |
| `revision/` | barthel-lab/KaryoScope-revisions | `b7d2092` | see `revision/README.md` |
| `benchmarks/a1-comparator-runtimes/` | barthel-lab/KaryoScope-revisions | `b7d2092` | |
| `benchmarks/tableS1/`, `benchmarks/figS6/` | cluster job directories (jobs 38434006, 38373939) | n/a | not previously in a repository |
| `benchmarks/figS6/plot_figS6.py` | barthel-lab/KaryoScope-paper | `main` | from `figures/fig3/data/`, with default paths made relative |
| `benchmarks/fig2/`, `benchmarks/edfig5/` | KaryoScope-archive working tree at `7c976e04`, and the Nature versions from KaryoScope-paper `main` | n/a | sha256 in `benchmarks/README.md` |
| `revision/a3-cross-genome/figure/` | barthel-lab/KaryoScope-paper | `main` | `scripts/make_build_figure.py`, `make_build_panels_ab.py` |
| `illustrator/` | barthel-lab/KaryoScope-paper | `main` | `scripts/nature_ai_export.py`, with the Arial Bold path set by `$ARIAL_BOLD` |
| `vendored/KaryoScope-ISCN/` (four figure files) | barthel-lab/KaryoScope-ISCN | `68a9816` | head of PR #87 (lowercase Fig. 3 panel letters) |
| `vendored/KaryoScope-ISCN/src/` | barthel-lab/KaryoScope-ISCN | `5507369` | 11 analysis modules, unmodified |
| `vendored/KaryoScope-conductor/` | barthel-lab/KaryoScope-conductor | `0732e96` | ten files |
| `vendored/KaryoScope-archive-resources/` | KaryoScope-archive pipeline | `7c976e04` | palettes, k-mer counts, and priority and hierarchy files |
| `vendored/KaryoScope-archive/bin_features.py` | barthel-lab/KaryoScope-archive (public) | `4434885` | binning script for the D4Z4 and centromere windows |

Commit `7c976e04` (the KaryoScope-archive pipeline commit named above) is not in
the public `barthel-lab/KaryoScope-archive` repository. It is reachable in the
private `barthel-lab/KaryoScope-archive-copy` and resolves through the public
`barthel-lab/KaryoScope-museum`.
