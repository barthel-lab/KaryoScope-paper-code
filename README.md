# KaryoScope-paper-code

Analysis and figure code for the manuscript *KaryoScope: rapid, alignment-free
sequence annotation for the pangenome era*. The KaryoScope tool itself is at
[barthel-lab/KaryoScope](https://github.com/barthel-lab/KaryoScope).

This repository is a snapshot: each directory is copied from the repository
and commit listed under Sources, with the changes stated in that directory's
README. It holds code and small result tables, not data.

## Figure and table map

Numbering as in the submitted manuscript.

| item | code | reproduced from this repository |
|---|---|---|
| Fig. 1b-e | `hg002/conf/hg002.yaml` -> `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py` | no: the published figure is an Illustrator composition, and its centromere file (`centromeres_v1coords.bed`) was not found |
| Fig. 2 | `benchmarks/fig2/figure2_benchmarking.R`; numbers: `revision/concordance/ks_concord.sh` | data yes; final artwork edited in Illustrator |
| Fig. 3 | `vendored/KaryoScope-ISCN/KaryoScope_assembly_contig_zoom_plot.py` via `vendored/KaryoScope-conductor`; configuration withheld (ROB coordinates, see below) | no: the annotations it was drawn from were not kept |
| Fig. 4 | `d4z4/make_figures.py` | yes, given the D4Z4 inputs (two are not public, see `d4z4/`) |
| Fig. 5a, 5c, 5d | `centromere/scripts/HPRC2_centromere_analysis.sh` | yes |
| Fig. 5b | none (FISH images) | n/a |
| ED Fig. 1-3 | none (Illustrator) | no |
| ED Fig. 4 | `revision/a3-cross-genome/figure/make_build_panels_ab.py` (a, b), `make_build_figure.py` (c); data `revision/a3-cross-genome/` | yes |
| ED Fig. 5 | `benchmarks/edfig5/supplementary_CHM13_benchmarking.nature.R` | yes |
| ED Fig. 6 | `benchmarks/figS6/` | yes |
| ED Fig. 7 | KaryoScope-ROBs figure scripts, withheld (ROB coordinates) | not in this snapshot |
| ED Fig. 8 | `d4z4/make_figures.py` | yes, as Fig. 4 |
| ED Fig. 9 | none (Illustrator); QC counts from `centromere/` steps 1-2 | no |
| ED Table 1 | `benchmarks/a1-comparator-runtimes/`, KaryoScope rows from `benchmarks/tableS1/` | scripts only; comparator logs were not kept |
| Supp. Fig. 1 | none (Illustrator) | no |
| Supp. Fig. 2 | `d4z4/make_figures.py` | yes, as Fig. 4 |
| Supp. Table 1 | `benchmarks/tableS1/` | yes, from the shipped run logs |
| Supp. Tables 2, 3 | `revision/a6-parameter-sweeps/`, `revision/a2-k-sensitivity/step4/` | sweep results included; the ROB rows are withheld |
| Supp. Notes 1, 2 | `revision/concordance/`, `revision/a4-smoothing-accounting/`, `benchmarks/` | scripts and result tables |

## Layout

```
README.md, LICENSE (GPL-3.0), CITATION.cff
environment/   one conda environment per analysis, exported from the environment it was run in
d4z4/          D4Z4 pipeline and figure code (Fig. 4, ED Fig. 8, Supp. Fig. 2)
centromere/    centromere clustering and gating (Fig. 5)
hg002/         HG002 assembly figure configuration (Fig. 1b-e)
benchmarks/    Supp. Table 1, ED Fig. 6, ED Table 1, Fig. 2 and ED Fig. 5 scripts
revision/      the revision analyses A2-A7 and the concordance harness
vendored/      the KaryoScope-ISCN and KaryoScope-conductor files the figure code imports
```

## Sources

| directory | repository | commit | notes |
|---|---|---|---|
| `d4z4/` | barthel-lab/HPRC-D4Z4 | `77e804b` | head of PR #1 (`figures/reproducible-layout`), not merged |
| `centromere/scripts/` | barthel-lab/KaryoScope-analysis | `90e2874` | branch `pangenome_structure_V2`, `scripts/` only |
| `hg002/` | barthel-lab/KaryoScope-HG002 | `27885ba` | `conf/`, `scripts/check_data.sh`, `requirements.txt` |
| `revision/` | barthel-lab/KaryoScope-revisions | `b7d2092` | see `revision/README.md`, Changes in this snapshot |
| `benchmarks/a1-comparator-runtimes/` | barthel-lab/KaryoScope-revisions | `b7d2092` | |
| `benchmarks/tableS1/`, `benchmarks/figS6/` | TGen cluster job directories `ks_slurm/tableS1_cold/`, `ks_slurm/threadsweep_cold/` | n/a | jobs 38434006 and 38373939; not previously in a repository |
| `benchmarks/figS6/plot_figS6.py` and its TSVs | barthel-lab/KaryoScope-paper | `6166acf` | `figures/fig3/data/` |
| `benchmarks/fig2/`, `benchmarks/edfig5/` | KaryoScope-archive working tree at `7c976e04` | n/a | uncommitted files; sha256 in `benchmarks/README.md` |
| `revision/a3-cross-genome/figure/` | barthel-lab/KaryoScope-paper | `bf6b6c3` (PR #35) | `scripts/make_build_figure.py`, `make_build_panels_ab.py`, with the Nature text sizes and lowercase panel letters |
| `vendored/KaryoScope-ISCN/` | barthel-lab/KaryoScope-ISCN | `4e42a98` | four files |
| `vendored/KaryoScope-conductor/` | barthel-lab/KaryoScope-conductor | `0732e96` | ten files |
| `vendored/KaryoScope-archive-resources/` | KaryoScope-archive pipeline | `7c976e04` | palettes, k-mer counts, priority and hierarchy files |

Commit `7c976e04` (the KaryoScope-archive pipeline commit named above) is not in the
public `barthel-lab/KaryoScope-archive` repository. It is reachable in the private
`barthel-lab/KaryoScope-archive-copy` and resolves through the public
`barthel-lab/KaryoScope-museum`.

Pending style-only updates (font sizes, panel-letter case, export format for
the Nature figure specification) are open as pull requests and are not in this
snapshot: HPRC-D4Z4 #2 (Fig. 4, ED Fig. 8, Supp. Fig. 2), KaryoScope-paper #35
(ED Fig. 6; the ED Fig. 4 scripts here are already the #35 versions),
KaryoScope-ISCN #87 (Fig. 3 lowercase panel letters) and KaryoScope-ROBs #2
(ED Fig. 7). The ED Fig. 5 styling is in `benchmarks/edfig5/` (`.nature.R`).
None of them changes plotted data.

## Withheld

The Robertsonian (ROB) cell-line assemblies are controlled-access (dbGaP
phs003920). Until the ruling on per-sample coordinates, everything that
holds ROB scaffold names or coordinates is left out: the KaryoScope-ROBs
configuration and figure scripts (Fig. 3 configuration, ED Fig. 7), A2 Step 3,
A6 `rob/`, and A7 Job 3 (details in `revision/README.md`). Their values come
from the dbGaP-approved assemblies. No sequence, read or per-sample BED is in
this repository.

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
