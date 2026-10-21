# KaryoScope revision analyses

Copied from `barthel-lab/KaryoScope-revisions` at commit `b7d2092`, with the
changes listed under "Changes in this snapshot" at the end of this file.

Analysis code and results generated in response to reviewer comments on the
KaryoScope manuscript. Each top-level directory is one analysis: the scripts that
produced it, the small result tables they emitted, and a `RESULTS.md` write-up
where one exists.

This repository is **code and results, not data**. Inputs are assemblies,
databases and annotation BEDs measured in tens of gigabytes; every script names
its inputs by path under the variables described in the repository README
(Paths), and `docs/INPUTS.md` records where they lived on the TGen cluster.

## Map: reviewer comment to analysis

| comment | question | analysis | headline |
|---|---|---|---|
| R1/R3 runtime | how does KaryoScope compare to the tools it is measured against? | [`../benchmarks/a1-comparator-runtimes`](../benchmarks/a1-comparator-runtimes) | RepeatMasker 43 h 54 m, HumAS-HMMER 53.10 h, Tiberius 2 h 16 m, against 13.71 min (KMC) / 7.91 min (HKS) for all six feature sets |
| R3-major-3 | is the annotation robust to *k*? | [`a2-k-sensitivity`](a2-k-sensitivity), [`step3`](a2-k-sensitivity/step3) (result table withheld, see below) | across *k* = 21–61 the Robertsonian breakpoint moves at most 900 bp; *k* = 11 is the control that breaks, and does |
| R3-major-3 | …and are the other biological applications? | [`a2-k-sensitivity/step4`](a2-k-sensitivity/step4) | **no.** The D4Z4 census fails at *k* = 51 (zero arrays) though it is salvageable by scaling the caller's thresholds with *k*; the centromere clustering keeps only 67.9% of calls between *k* = 21 and 51. Novel fraction and the divergence series are measured here too |
| R3-major-1 | does this work on genomes other than human? | [`a3-cross-genome`](a3-cross-genome) | Arabidopsis Col-CEN and Ey15-2; label composition per feature set |
| R2/R3 smoothing | what does smoothing actually do, and to how much? | [`a4-smoothing-accounting`](a4-smoothing-accounting) | novel and multilabel accounting per feature set; 2.5643% of positions carry `novel` before smoothing, identical across all six sets |
| R3 divergence | how far from the reference does it still work? | [`a5-divergence-simulation`](a5-divergence-simulation) | 13 mutated-genome conditions, 100% coverage |
| Minor 1 | a parameter table with justifications, and a sensitivity note | [`a6-parameter-sweeps`](a6-parameter-sweeps) | five sweeps; **no published number moves** |
| 2026-09-07 review | four factual checks and one mechanism check before resubmission | [`a7-review-checks`](a7-review-checks) | Ey15-2 is **18.66 % unrepresented** before smoothing (the 105 bp was the imputation residue); the gene-set percentages in the Results are **label-permuted** (exon is 11.5 %, not 49.4 %); smoothing is **exactly** reverse-complement symmetric on both backends; the acrocentric flank rule holds once "flank" means the nearest same-branch label; D4Z4 same-contig is by construction |
| — | the concordance harness the above rely on | [`concordance`](concordance) | reproduces the manuscript exactly on both HG002 haplotypes |

## The parameter sweeps, in one table

Everything in [`a6-parameter-sweeps`](a6-parameter-sweeps). Each has its own
`RESULTS.md` with the full numbers and a verdict for the supplementary table.

| sweep | parameter(s) | result |
|---|---|---|
| [`maxgap`](a6-parameter-sweeps/maxgap) | smoothing `max_gap` = 1000 bp | concordance is **flat to 0.0136 points** from 0 to 10,000 bp, so concordance cannot be what selected the value; justified instead from the gap-length distribution (99.4% of novel bases in runs <= 1 kb, median novel run exactly 31 bp = *k*) |
| [`tar1-its`](a6-parameter-sweeps/tar1-its) | TAR1/ITS merge distance = 1 kb | **no knee** in fragment count; the evidence is the bridging ceiling, stable to d = 2000 then doubling, so 1 kb sits in a wide safe band |
| [`rob`](a6-parameter-sweeps/rob) (result table withheld, see below) | ROB window / coverage / dominance | the three parameters named in the response letter control the ISCN *interval*, not the breakpoint: **the published breakpoint is identical across all 39 configurations**. The two that do reach it move it by at most 9.1 kb over a 100-fold range |
| [`centromere-clustering`](a6-parameter-sweeps/centromere-clustering) | `--sil-threshold` 0.5, `--centroid-sd` 5 | **13 chromosomes retain a split at all 30 settings**; the four chromosomes in the allele heatmap are invariant at all 30; 0.5 sits on a plateau with 0.6 |
| [`d4z4`](a6-parameter-sweeps/d4z4) | ten merge thresholds | **the published repeat count is 28,580 in all 51 runs** — bit-identical under every threshold; five of the ten are completely inert, and the one exposed parameter is undisclosed |

## Conventions used throughout

These are worth knowing before reading or re-running anything.

* **Runtimes are attributed to the node class.** Everything quoted is one Intel
  Xeon Gold 6342 node (Ice Lake, 2.80 GHz), matching Table S1. Jobs pin it with
  `--constraint=2.80GHz&Ice_Lake`.
* **Coordinates are stated in one frame.** Scaffolded and original-assembly
  coordinates differ on any flipped contig, and mixing them produces a
  plausible-looking megabase-scale error rather than a crash. The A2 write-ups
  name the frame in every table.
* **Never glob across the presmoothed/smoothed pair.** `*.smoothed.bed.gz` also
  matches `presmoothed`, and it sorts *first*, so `ls … | head -1` silently
  selects the wrong input. Scripts here name files exactly.
* **A control that must reproduce a known answer runs in every sweep.** The *k*
  sweep runs *k* = 31 through the same driver and requires it to reproduce the
  shipped pipeline to the base; the centromere sweep runs the published (0.5, 5)
  through all of steps 4–7. A sweep whose control fails is discarded, not
  interpreted.

## External code

Three analyses drive code that lives in other repositories. It is **not vendored
here** — it is referenced by repository and commit so that the version used is
recorded without duplicating it:

| analysis | code | version used |
|---|---|---|
| `a6-parameter-sweeps/centromere-clustering` | `barthel-lab/KaryoScope-analysis`, branch `pangenome_structure_V2`, `scripts/HPRC2_centromere_analysis.sh` and the five scripts it calls | `90e2874`, run unmodified |
| `a6-parameter-sweeps/d4z4` | `barthel-lab/HPRC-D4Z4`, `analyze_d4z4_arrays.py` | one line changed: `base_dir` read from `D4Z4_BASE_DIR`. The copy in `../d4z4/` already contains this change |
| `a1-comparator-runtimes` | RepeatMasker 4.1.6 / Dfam 3.8, HumAS-HMMER for AnVIL, Tiberius | the HumAS-HMMER `--dna` fix is recorded in that directory's README |

## Layout

```
a2-k-sensitivity/           the k sweep (Steps 0-2)
a2-k-sensitivity/step3/     the Robertsonian breakpoint vs k (ROB configuration needed)
a2-k-sensitivity/step4/     novel fraction, divergence, D4Z4 and centromere clustering vs k
a3-cross-genome/            Arabidopsis Col-CEN and Ey15-2
a4-smoothing-accounting/    novel and multilabel accounting per feature set
a5-divergence-simulation/   mutated-genome conditions
a6-parameter-sweeps/        the five sensitivity sweeps behind the parameter table
a7-review-checks/           the five checks raised by the 2026-09-07 independent review
concordance/                the harness the accuracy numbers rest on
docs/INPUTS.md              where each input lived on the TGen cluster
```

Each analysis directory holds `scripts/` (as run, including SLURM directives) and
`results/` (the small tables). Where a directory has a `RESULTS.md`, read that
first.

## Changes in this snapshot

- `a1-comparator-runtimes/` is in `../benchmarks/`.
- Not copied: `manuscript-tables/` (superseded copies of Supplementary Tables 2
  and 3), `figures/figS6/` (the August thread sweep, superseded by
  `../benchmarks/figS6/`), `docs/briefs/` and `docs/write-ups/` (internal task
  briefs and hand-offs).
- **The Robertsonian (ROB) analyses read their per-sample values from the ROB
  configuration.** The scripts of A2 Step 3 (`a2-k-sensitivity/step3/`), the A6
  ROB sweep (`a6-parameter-sweeps/rob/`) and A7 Job 3 (`a7-review-checks/job3/`)
  read every ROB scaffold name, coordinate and sample-specific value from
  `$ROBS_CONFIG`, sections `a2_step3`, `a6_rob` and `a7_job3` of
  `../robs/conf/config.template.yaml`; at `b7d2092` these values are written in
  the scripts. With the filled configuration (dbGaP phs003920) each script gives
  byte-identical output to the original on the same inputs. `a2-k-sensitivity/step3/README.md`
  replaces the Step 3 part of `a2-k-sensitivity/README.md`, and
  `a6-parameter-sweeps/rob/README.md` summarises `rob/RESULTS.md` without
  coordinates.
- **Withheld pending the dbGaP ruling on the ROB assemblies:** the result
  tables and write-ups that hold ROB coordinates or scaffold names:
  `a2-k-sensitivity/results/sst1_breakpoints_by_k.tsv`, the Step 3 part of
  `a2-k-sensitivity/README.md`, `a6-parameter-sweeps/rob/results/rob_param_sweep.tsv`
  and `rob/RESULTS.md`, the A7 Job 3 provenance logs, and the ROB rows of
  `a7-review-checks/results/rc_symmetry_*.tsv` (the HG002 rows are complete).
  The Job 3 paragraph in `a7-review-checks/RESULTS.md` does not name the fusion
  scaffold.
- `a6-parameter-sweeps/d4z4/base_dir.patch` is dropped (see External code).
- Cluster paths are replaced by the variables in the repository README (Paths).
  Result tables that recorded input paths (`a2-k-sensitivity/step4/results/*_windows.tsv`)
  carry the same substitution.
