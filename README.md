# KaryoScope-paper-code

Analysis and figure code for the manuscript *Rapid alignment-free annotation for
the pangenome era with KaryoScope*. The bioRxiv preprint carries the original
title, *KaryoScope: rapid, alignment-free sequence annotation for the pangenome
era*. The KaryoScope tool itself is at
[barthel-lab/KaryoScope](https://github.com/barthel-lab/KaryoScope).

This repository is a snapshot: each directory is copied from the repository
and commit listed under Sources, with the changes stated in that directory's
README. It holds code and small result tables, not data.

## Figure and table map

Numbering as in the revised manuscript.

| item | code | reproduced from this repository |
|---|---|---|
| Fig. 1b-e | `hg002/conf/hg002.yaml` -> `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py`; centromere windows `hg002/scripts/make_centromeres_v1coords.py`; submitted file from the Illustrator composition with `illustrator/nature_ai_export.py --relabel C=E,E=C` | panels yes; the layout is an Illustrator composition. Its centromere window file (`centromeres_v1coords.bed`) is rebuilt from the v1 and v2 centromere BEDs (`hg002/data/`); all 46 window lengths match Fig. 1b |
| Fig. 2 | `benchmarks/fig2/figure2_benchmarking.nature.R` (submitted figure; the original `figure2_benchmarking.R` beside it); numbers: `revision/concordance/ks_concord.sh` | yes. The submitted PDF was drawn on macOS with the Quartz device; elsewhere the script falls back to `cairo_pdf` (same text and data, different glyph spacing) |
| Fig. 3 | `robs/` (configuration template) -> `vendored/KaryoScope-ISCN/KaryoScope_assembly_contig_zoom_plot.py` via `vendored/KaryoScope-conductor` | with dbGaP access, from fresh v1.1.0 annotations (see `robs/README.md`): on the tested sample the render is byte-identical to the submitted figure before a 4-unit headroom edit, which `robs/scripts/fig3_headroom.py` applies |
| Fig. 4 | `d4z4/make_figures.py` | yes, given the D4Z4 inputs (two are not public, see `d4z4/`) |
| Fig. 5a, 5c, 5d | `centromere/scripts/HPRC2_centromere_analysis.sh` | yes |
| Fig. 5b | none (FISH images); the submitted Fig. 5 file from the Illustrator composition with `illustrator/nature_ai_export.py` | n/a |
| ED Fig. 1-3 | Illustrator; submitted files from `illustrator/nature_ai_export.py` (ED Fig. 2 with `--max-pt 7`) | artwork no |
| ED Fig. 4 | `revision/a3-cross-genome/figure/make_build_panels_ab.py` (a, b), `make_build_figure.py` (c); data `revision/a3-cross-genome/` | yes |
| ED Fig. 5 | `benchmarks/edfig5/supplementary_CHM13_benchmarking.nature.R` | yes |
| ED Fig. 6 | `benchmarks/figS6/` | yes |
| ED Fig. 7 | `robs/scripts/` (configuration template) | with dbGaP access, from fresh v1.1.0 annotations (see `robs/README.md`) |
| ED Fig. 8 | `d4z4/make_figures.py` | yes, as Fig. 4 |
| ED Fig. 9 | Illustrator, exported with `illustrator/nature_ai_export.py --max-pt 7 --anchor center --drop-small`; QC table `centromere/qc/`, counts from `centromere/` steps 1-2 | artwork no; the QC table and the 6,138 / 4,650 counts yes, from the public HPRC QC calls |
| ED Table 1 | `benchmarks/a1-comparator-runtimes/`, KaryoScope rows from `benchmarks/tableS1/` | scripts only; comparator logs were not kept |
| Supp. Fig. 1 | Illustrator, exported with `illustrator/nature_ai_export.py` | artwork no |
| Supp. Fig. 2 | `d4z4/make_figures.py` | yes, as Fig. 4 |
| Supp. Table 1 | `benchmarks/tableS1/` | yes, from the shipped run logs |
| Supp. Tables 2, 3 | `revision/a6-parameter-sweeps/` (ROB block of Supp. Table 2: `rob/`), `revision/a2-k-sensitivity/step4/` | sweep results included; the ROB rows with dbGaP access, from the ROB configuration (the ROB result table is withheld) |
| Supp. Notes 1, 2 | `revision/concordance/`, `revision/a4-smoothing-accounting/`, `benchmarks/` | scripts and result tables |

## Versions

- **KaryoScope v1.1.0** (tag `4c09ce5`, doi:10.5281/zenodo.20657817): the
  annotations behind Figs. 1-2 (HG002 v1.1, `KS_human_CHM13_v2`) and the D4Z4
  and centromere inputs (HPRC Release 2, `KS_human_CHM13`, presmoothed). These
  were first made with the KaryoScope-archive Snakemake pipeline (`7c976e04`
  for HG002, `0293e85` for HPRC); v1.1.0 reproduces them byte for byte.
- **Binning of the D4Z4 and centromere windows** (100 bp and 100 kb):
  `vendored/KaryoScope-archive/bin_features.py` (KaryoScope-archive `4434885`).
- **KaryoScope v2.3.0** (tag `124f84a`, doi:10.5281/zenodo.21852728) with
  **hks 0.4.0** (HKS `04bbc52`): the benchmarks (Supp. Table 1, ED Fig. 6), the
  tool comparison (ED Table 1), the k-length (A2), divergence (A5) and
  cross-genome (A3) analyses. Supp. Table 1 and ED Fig. 6 ran on 2026-10-12
  and 2026-10-13 from a clean worktree of the v2.3.0 tag. The other analyses ran
  between 2026-08-11 and 2026-09-03 from a checkout seven or eight commits past
  the tag (`052055a`, `501d5fd`; it reports "2.3.0"). Those commits change only
  the `karyotype`, `scaffold` and haplotype-label code, not `annotate`, `bin` or
  `build`. One exception: the A2 databases at k = 11-61 were built on
  2026-08-03, before the v2.3.0 tag, with `karyoscope build` from `9ccba0b`
  (v2.1.0 + 3 commits) and hks `c5fecf8` (the 0.4.0 version before the upstream
  merge); the annotations at each k used the v2.3.0 code with hks `04bbc52`.
- **Robertsonian figures** (Fig. 3, ED Fig. 7): the assemblies were first
  annotated, scaffolded and binned with the KaryoScope-archive pipeline (commit
  not recorded). KaryoScope v1.1.0 (annotate and scaffold) with `karyoscope bin`
  (v1.1.0 or v2.3.0) reproduces the submitted Fig. 3 byte for byte on the tested
  sample; see `robs/README.md`.

## Layout

```
README.md, LICENSE (GPL-3.0), CITATION.cff
environment/   one conda environment per analysis, exported from the environment it was run in
d4z4/          D4Z4 pipeline and figure code (Fig. 4, ED Fig. 8, Supp. Fig. 2)
centromere/    centromere clustering and gating (Fig. 5); qc/: the assembly QC table (ED Fig. 9)
hg002/         HG002 assembly figure configuration (Fig. 1b-e)
robs/          Robertsonian translocation figure code and configuration template (Fig. 3, ED Fig. 7)
benchmarks/    Supp. Table 1, ED Fig. 6, ED Table 1, Fig. 2 and ED Fig. 5 scripts
illustrator/   export of the Illustrator-composed figures to the Nature specification
revision/      the revision analyses A2-A7 and the concordance harness
vendored/      the KaryoScope-ISCN and KaryoScope-conductor files the figure code imports, and the
               KaryoScope-archive binning script (bin_features.py)
```

## Sources

| directory | repository | commit | notes |
|---|---|---|---|
| `d4z4/` | barthel-lab/HPRC-D4Z4 | `77e804b`; figure code `52ed8cc` | heads of PR #1 (`figures/reproducible-layout`) and PR #2 (`nature/ed-spec`, the Nature figure size), not merged |
| `robs/` | barthel-lab/KaryoScope-ROBs | `1188850` | head of PR #2 (`nature/ed-spec`, on `fig5/assemble-panels`), not merged; configuration replaced by a template |
| `centromere/scripts/` | barthel-lab/KaryoScope-analysis | `90e2874` | branch `pangenome_structure_V2`, `scripts/` only |
| `hg002/` | barthel-lab/KaryoScope-HG002 | `27885ba` | `conf/`, `scripts/check_data.sh`, `requirements.txt` |
| `hg002/scripts/make_centromeres_v1coords.py`, `hg002/data/` | written for this snapshot | n/a | rebuilds the Fig. 1b centromere window file |
| `centromere/qc/` | KaryoScope-archive `0bf3a7d` (`workflow/rules/pangenome.smk`, `get_regions.py`) and two files that were never committed (`aggregate_master_qc.py`, `aggregate_qc_v2.awk`) | n/a | originals in `centromere/qc/original/`; the runnable scripts were written for this snapshot and reproduce the 2025 tables byte for byte |
| `revision/` | barthel-lab/KaryoScope-revisions | `b7d2092` | see `revision/README.md`, Changes in this snapshot |
| `benchmarks/a1-comparator-runtimes/` | barthel-lab/KaryoScope-revisions | `b7d2092` | |
| `benchmarks/tableS1/`, `benchmarks/figS6/` | TGen cluster job directories `ks_slurm/tableS1_cold/`, `ks_slurm/threadsweep_cold/` | n/a | jobs 38434006 and 38373939; not previously in a repository |
| `benchmarks/figS6/plot_figS6.py` and its TSVs | barthel-lab/KaryoScope-paper | `605aa7f` (`main`) | `figures/fig3/data/`; the script is the Nature-size version (PR #35) |
| `benchmarks/fig2/`, `benchmarks/edfig5/` | KaryoScope-archive working tree at `7c976e04` | n/a | uncommitted files; sha256 in `benchmarks/README.md` |
| `revision/a3-cross-genome/figure/` | barthel-lab/KaryoScope-paper | `605aa7f` (`main`) | `scripts/make_build_figure.py`, `make_build_panels_ab.py`, with the Nature text sizes and lowercase panel letters |
| `vendored/KaryoScope-ISCN/` | barthel-lab/KaryoScope-ISCN | `68a9816` | head of PR #87 (lowercase Fig. 3 panel letters), not merged; four files |
| `vendored/KaryoScope-conductor/` | barthel-lab/KaryoScope-conductor | `0732e96` | ten files |
| `vendored/KaryoScope-archive-resources/` | KaryoScope-archive pipeline | `7c976e04` | palettes, k-mer counts, priority and hierarchy files |
| `vendored/KaryoScope-archive/bin_features.py` | barthel-lab/KaryoScope-archive (public) | `4434885` | the binning script for the D4Z4 and centromere windows |

Commit `7c976e04` (the KaryoScope-archive pipeline commit named above) is not in the
public `barthel-lab/KaryoScope-archive` repository. It is reachable in the private
`barthel-lab/KaryoScope-archive-copy` and resolves through the public
`barthel-lab/KaryoScope-museum`.

## Withheld

The Robertsonian (ROB) cell-line assemblies are controlled-access (dbGaP
phs003920). ROB scaffold names and coordinates are not in this repository.
The code that needs them reads them from the ROB configuration
(`$ROBS_CONFIG`): `robs/conf/config.template.yaml` has a placeholder for each
value and `robs/README.md` gives the steps to fill it. This covers the ROB
figures (`robs/`), A2 Step 3 (`revision/a2-k-sensitivity/step3/`), the A6 ROB
sweep (`revision/a6-parameter-sweeps/rob/`) and A7 Job 3
(`revision/a7-review-checks/job3/`). Withheld are only the result tables and
write-ups that hold ROB values (listed in `revision/README.md`, Changes in this
snapshot); a run with dbGaP access regenerates them. No sequence, read or
ROB per-sample BED is in this repository.

## Paths

Scripts that ran on the TGen cluster name their inputs under these variables
(set them before running; a shell script stops with a message if one is unset,
and a Python script falls back to the literal text `$NAME`, so the first file
access fails with that path):

| variable | what it held on the cluster |
|---|---|
| `KS_WORK` | the working area (analysis outputs, staged inputs, job directories) |
| `KS_PROJECTS` | project data (KaryoScope results for HG002, CHM13 and the HPRC pangenome) |
| `KS_REFERENCES` | reference assemblies (CHM13, HG002 v1.1, HPRC Release 2) |
| `KS_SOFTWARE` | software (conda environments, the KaryoScope checkout and its database registry) |
| `KS_HOME` | the home directory (R and tool binaries) |
| `KS_GUARRACINO` | the directory holding the ROB assemblies |
| `KS_RESOURCES` | `vendored/KaryoScope-archive-resources/resources` (figure palettes and trees) |
| `D4Z4_BASE_DIR` | the D4Z4 working directory (`analyze_d4z4_arrays.py`; default `.`) |

`#SBATCH --output` lines write to the submission directory.

## Inputs

- HPRC Release 2 assemblies: the HPRC data release (public).
- HG002 v1.1 and T2T-CHM13v2.0: public.
- KaryoScope databases (`KS_human_CHM13`, `KS_human_CHM13_v2`, `HKS_*`):
  installed with the KaryoScope tool (see its README).
- Intermediate annotations (BEDs, tens of GB) are regenerated with KaryoScope;
  each directory's README names the files it reads.

## Licence

GPL-3.0-or-later (`LICENSE`).
