# KaryoScope-paper-code

Code to reproduce the figures and analyses in *Rapid alignment-free annotation
for the pangenome era with KaryoScope* (bioRxiv preprint title: *KaryoScope:
rapid, alignment-free sequence annotation for the pangenome era*).

KaryoScope itself is a separate tool:
[barthel-lab/KaryoScope](https://github.com/barthel-lab/KaryoScope). This
repository contains the analysis and plotting code built around it, and small
result tables, but no sequence data.

## Before you start

1. **Install KaryoScope and its databases** following the
   [KaryoScope README](https://github.com/barthel-lab/KaryoScope). The analyses
   use v1.1.0 and v2.3.0 (see [Versions](#versions)).
2. **Create the conda environment** for the analysis you want to run, from
   `environment/` (one file per analysis: `figures.yml`, `centromere.yml`,
   `d4z4.yml`, and `benchmarks-r.yml`).
3. **Get the input data.** The inputs are public (T2T-CHM13v2.0, HG002 v1.1,
   HPRC Release 2, and its QC calls), except the three Robertsonian
   translocation cell-line assemblies, which are controlled-access (see
   [Controlled-access data](#controlled-access-data)), and two D4Z4 inputs (see
   `d4z4/README.md`).
4. **Set the path variables** the scripts read (see [Paths](#paths)).

Each directory has a README with its own run steps and inputs.

## Where each figure and table comes from

Numbering follows the revised manuscript. ED is Extended Data, and Supp. is
Supplementary.

| item | directory | reproducible from here |
|---|---|---|
| Fig. 1a | — | no: schematic |
| Fig. 1b–e | `hg002/` | panels yes (final layout assembled in Illustrator¹) |
| Fig. 2 | `benchmarks/fig2/` | yes |
| Fig. 3 | `robs/` | yes, with dbGaP access² |
| Fig. 4 | `d4z4/` | yes³ |
| Fig. 5a, c, d | `centromere/` | yes |
| Fig. 5b | — | no: FISH images |
| ED Figs. 1–3, Supp. Fig. 1 | — | no: Illustrator figures¹ |
| ED Fig. 4 | `revision/a3-cross-genome/` | yes |
| ED Fig. 5 | `benchmarks/edfig5/` | yes |
| ED Fig. 6 | `benchmarks/figS6/` | yes |
| ED Fig. 7 | `robs/` | yes, with dbGaP access² |
| ED Fig. 8, Supp. Fig. 2 | `d4z4/` | yes³ |
| ED Fig. 9 | `centromere/qc/` | QC table and counts yes (layout assembled in Illustrator¹) |
| ED Table 1 | `benchmarks/a1-comparator-runtimes/`, `benchmarks/tableS1/` | scripts only (the comparator run logs were not kept) |
| Supp. Table 1 | `benchmarks/tableS1/` | yes, from the included run logs |
| Supp. Tables 2, 3 | `revision/a6-parameter-sweeps/`, `revision/a2-k-sensitivity/` | yes (the Robertsonian rows need dbGaP access²) |
| Supp. Notes 1, 2 | `revision/concordance/`, `revision/a4-smoothing-accounting/`, `benchmarks/` | yes |

¹ The figures composed in Illustrator were exported to the journal format with
`illustrator/nature_ai_export.py`. The options used for each figure are listed
in `illustrator/README.md`. The `.ai` files are not in this repository.
² See [Controlled-access data](#controlled-access-data) and `robs/README.md`.
³ Two D4Z4 inputs (the HMM repeat counts and the pLAM BLAST table) are not
public. See `d4z4/README.md`.

The analyses added during revision (A1–A7: runtime comparison, *k*-length,
cross-genome, smoothing, divergence, parameter sweeps, and review checks) are
indexed in `revision/README.md`.

## Versions

| what | KaryoScope version |
|---|---|
| Annotations behind Figs. 1–3, ED Fig. 7, and the D4Z4 and centromere inputs | v1.1.0 ([doi:10.5281/zenodo.20657817](https://doi.org/10.5281/zenodo.20657817)) |
| Runtime and memory benchmarks, tool comparison, and the *k*-length, divergence, and cross-genome analyses | v2.3.0 ([doi:10.5281/zenodo.21852728](https://doi.org/10.5281/zenodo.21852728)) with hks 0.4.0 |

The D4Z4 and centromere windows were binned with
`vendored/KaryoScope-archive/bin_features.py`, and the Robertsonian windows with
`karyoscope bin` (v1.1.0). Exact commits and run dates are in
[PROVENANCE.md](PROVENANCE.md), with two exceptions noted there: the *k*-length
databases were built before v2.3.0, and the annotations were first made with an
earlier pipeline that v1.1.0 reproduces.

## Repository layout

```
hg002/         Fig. 1 (HG002 assembly figure)
benchmarks/    Fig. 2, ED Figs. 5-6, ED Table 1, Supp. Table 1
robs/          Fig. 3, ED Fig. 7 (Robertsonian translocations)
d4z4/          Fig. 4, ED Fig. 8, Supp. Fig. 2
centromere/    Fig. 5 (clustering) and qc/ (assembly QC table, ED Fig. 9)
revision/      analyses A1-A7 added during revision, and the concordance harness
illustrator/   export of the figures composed in Illustrator
vendored/      code copied from other repositories that the scripts import
environment/   conda environments, one per analysis
```

## Controlled-access data

The three Robertsonian translocation cell-line assemblies (GM03417, GM03786, and
GM04890) are in dbGaP (phs003920.v1.p1). This repository contains no sequence,
annotation, scaffold name, or coordinate from them.

Code that needs those values reads them from a configuration file that you fill
in: `robs/conf/config.template.yaml`, passed as `$ROBS_CONFIG`. The steps from
dbGaP access to a finished figure are in `robs/README.md`. The same
configuration is used for the Robertsonian parts of the revision analyses
(`revision/a2-k-sensitivity/step3/`, `revision/a6-parameter-sweeps/rob/`, and
`revision/a7-review-checks/job3/`). Their result tables are not included, and
can be regenerated with dbGaP access.

## Paths

The scripts were run on a compute cluster and locate their inputs through these
variables. Shell scripts stop with a message if one is unset.

| variable | points to |
|---|---|
| `KS_WORK` | a working directory for outputs and staged inputs |
| `KS_PROJECTS` | KaryoScope results for HG002, CHM13, and the HPRC pangenome |
| `KS_REFERENCES` | reference assemblies (CHM13, HG002 v1.1, and HPRC Release 2) |
| `KS_SOFTWARE` | conda environments, the KaryoScope checkout, and its database registry |
| `KS_HOME` | home directory (R and tool binaries) |
| `KS_GUARRACINO` | the Robertsonian assemblies (dbGaP) |
| `KS_RESOURCES` | `vendored/KaryoScope-archive-resources/resources` |
| `KS_ISCN_SRC` | `vendored/KaryoScope-ISCN/src` |
| `D4Z4_BASE_DIR` | the D4Z4 working directory (default `.`) |

## Citation and licence

Cite the manuscript (see `CITATION.cff`). The code is licensed GPL-3.0-or-later
(`LICENSE`).
