# KaryoScope-paper-code

Code to reproduce the figures and analyses in *Rapid alignment-free annotation
for the pangenome era with KaryoScope* (bioRxiv preprint title: *KaryoScope:
rapid, alignment-free sequence annotation for the pangenome era*).

KaryoScope itself is a separate tool:
[barthel-lab/KaryoScope](https://github.com/barthel-lab/KaryoScope). This
repository holds the analysis and plotting code built around it, plus small
result tables. It holds no sequence data.

## Before you start

1. **Install KaryoScope and its databases** following the
   [KaryoScope README](https://github.com/barthel-lab/KaryoScope). The analyses
   use v1.1.0 and v2.3.0 (see [Versions](#versions)).
2. **Create the conda environment** for the analysis you want to run, from
   `environment/` (one file per analysis: `figures.yml`, `centromere.yml`,
   `d4z4.yml`, `benchmarks-r.yml`).
3. **Get the input data.** Almost everything is public (T2T-CHM13v2.0, HG002
   v1.1, HPRC Release 2 and its QC calls). The three Robertsonian translocation
   cell lines are controlled-access; see [Controlled-access data](#controlled-access-data).
4. **Set the path variables** the scripts read (see [Paths](#paths)).

Each directory has a README with its own run steps and inputs.

## Where each figure and table comes from

Numbering follows the revised manuscript.

| item | directory | reproducible from here |
|---|---|---|
| Fig. 1b–e | `hg002/` | panels yes; final layout assembled in Illustrator¹ |
| Fig. 2 | `benchmarks/fig2/` | yes |
| Fig. 3 | `robs/` | yes, with dbGaP access² |
| Fig. 4 | `d4z4/` | yes³ |
| Fig. 5a, c, d | `centromere/` | yes |
| Fig. 5b | — | no: FISH images |
| ED Fig. 1–3, Supp. Fig. 1 | — | no: Illustrator figures¹ |
| ED Fig. 4 | `revision/a3-cross-genome/` | yes |
| ED Fig. 5 | `benchmarks/edfig5/` | yes |
| ED Fig. 6 | `benchmarks/figS6/` | yes |
| ED Fig. 7 | `robs/` | yes, with dbGaP access² |
| ED Fig. 8, Supp. Fig. 2 | `d4z4/` | yes³ |
| ED Fig. 9 | `centromere/qc/` | QC table and counts yes; layout assembled in Illustrator¹ |
| ED Table 1 | `benchmarks/a1-comparator-runtimes/`, `benchmarks/tableS1/` | scripts only (comparator run logs were not kept) |
| Supp. Table 1 | `benchmarks/tableS1/` | yes, from the included run logs |
| Supp. Tables 2, 3 | `revision/a6-parameter-sweeps/`, `revision/a2-k-sensitivity/` | yes; the Robertsonian rows need dbGaP access² |
| Supp. Notes 1, 2 | `revision/concordance/`, `revision/a4-smoothing-accounting/`, `benchmarks/` | yes |

¹ The Illustrator-composed figures were exported to the journal format with
`illustrator/nature_ai_export.py`; its README lists the options used per figure.
The `.ai` files are not in this repository.
² See [Controlled-access data](#controlled-access-data) and `robs/README.md`.
³ Two D4Z4 inputs (the HMM repeat counts and the pLAM BLAST table) are not
public; see `d4z4/README.md`.

The revision analyses (A1–A7: runtime comparison, *k*-length, cross-genome,
smoothing, divergence, parameter sweeps, review checks) are indexed in
`revision/README.md`.

## Versions

| what | KaryoScope version |
|---|---|
| Annotations behind Figs. 1–3, ED Fig. 7, and the D4Z4 and centromere inputs | v1.1.0 ([doi:10.5281/zenodo.20657817](https://doi.org/10.5281/zenodo.20657817)) |
| Benchmarks, tool comparison, *k*-length, divergence and cross-genome analyses | v2.3.0 ([doi:10.5281/zenodo.21852728](https://doi.org/10.5281/zenodo.21852728)) with hks 0.4.0 |

The D4Z4 and centromere windows were binned with
`vendored/KaryoScope-archive/bin_features.py`. Exact commits, run dates and the
two exceptions (the *k*-length databases and the original pipeline runs) are in
[PROVENANCE.md](PROVENANCE.md).

## Repository layout

```
hg002/         Fig. 1 (HG002 assembly figure)
benchmarks/    Fig. 2, ED Figs. 5-6, ED Table 1, Supp. Table 1
robs/          Fig. 3, ED Fig. 7 (Robertsonian translocations)
d4z4/          Fig. 4, ED Fig. 8, Supp. Fig. 2
centromere/    Fig. 5 (clustering) and qc/ (assembly QC table, ED Fig. 9)
revision/      revision analyses A1-A7 and the concordance harness
illustrator/   export of the Illustrator-composed figures
vendored/      code copied from other repositories that the scripts import
environment/   conda environments, one per analysis
```

## Controlled-access data

The three Robertsonian translocation cell-line assemblies (GM03417, GM03786,
GM04890) are in dbGaP (phs003920). This repository contains no sequence,
annotation, scaffold name or coordinate from them.

Code that needs those values reads them from a configuration file you fill in
yourself: `robs/conf/config.template.yaml`, passed as `$ROBS_CONFIG`.
`robs/README.md` gives the steps, from dbGaP access to a finished figure. The
same configuration drives the Robertsonian parts of the revision analyses
(`revision/a2-k-sensitivity/step3/`, `revision/a6-parameter-sweeps/rob/`,
`revision/a7-review-checks/job3/`). Their result tables are not included; a run
with dbGaP access regenerates them.

## Paths

The scripts ran on a compute cluster and find their inputs through these
variables. Shell scripts stop with a message if one is unset.

| variable | points to |
|---|---|
| `KS_WORK` | a working directory for outputs and staged inputs |
| `KS_PROJECTS` | KaryoScope results for HG002, CHM13 and the HPRC pangenome |
| `KS_REFERENCES` | reference assemblies (CHM13, HG002 v1.1, HPRC Release 2) |
| `KS_SOFTWARE` | software: conda environments, the KaryoScope checkout and its database registry |
| `KS_HOME` | home directory (R and tool binaries) |
| `KS_GUARRACINO` | the Robertsonian assemblies |
| `KS_RESOURCES` | `vendored/KaryoScope-archive-resources/resources` |
| `KS_ISCN_SRC` | `vendored/KaryoScope-ISCN/src` |
| `D4Z4_BASE_DIR` | the D4Z4 working directory (default `.`) |

## Citation and licence

Cite the manuscript (see `CITATION.cff`). Code is GPL-3.0-or-later (`LICENSE`).
