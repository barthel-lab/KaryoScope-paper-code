# KaryoScope revision analyses

Copied from `barthel-lab/KaryoScope-revisions` at commit `b7d2092`, with the
changes listed under "Changes in this snapshot" at the end of this file.

Analysis code and results generated in response to reviewer comments on the
KaryoScope manuscript. Each top-level directory is one analysis: the scripts that
produced it, the small result tables they wrote, and a `RESULTS.md` write-up
where one exists.

This directory contains **code and results, not data**. The inputs are
assemblies, databases, and annotation BEDs of tens of gigabytes. Every script
names its inputs by path under the variables described in the repository README
(Paths), and where they were on the TGen cluster is recorded in
`docs/INPUTS.md`.

## Map: reviewer comment to analysis

R1, R2, and R3 are reviewers 1, 2, and 3.

| comment | question | analysis | headline |
|---|---|---|---|
| R1/R3 runtime | how does the runtime of KaryoScope compare with that of the specialized tools? | [`../benchmarks/a1-comparator-runtimes`](../benchmarks/a1-comparator-runtimes) | RepeatMasker 43 h 54 m, HumAS-HMMER 53.10 h, Tiberius 2 h 16 m, compared with 11.5 min (KMC) and 6.3 min (HKS) for all six feature sets (*n* = 3 means, Supplementary Table 1) |
| R3-major-3 | is the annotation robust to *k*? | [`a2-k-sensitivity`](a2-k-sensitivity), [`step3`](a2-k-sensitivity/step3) (result table withheld, see below) | across *k* = 21–61 the Robertsonian breakpoint moves at most 900 bp, and *k* = 11 is the control expected to fail, which it does |
| R3-major-3 | …and are the other biological applications? | [`a2-k-sensitivity/step4`](a2-k-sensitivity/step4) | **not without changes.** With its thresholds unchanged, the D4Z4 array caller finds no arrays at *k* = 51 and half as many at *k* = 21. With its two size thresholds scaled to *k*, it calls 630 arrays at *k* = 21 and 51 (629 at *k* = 31). The major or minor calls at *k* = 31 are reproduced for 76.4% (*k* = 21) and 85.9% (*k* = 51) of centromere haplotypes; 67.9% keep their call between *k* = 21 and 51, and 2 of the 4 FISH-validated subtypes are lost away from *k* = 31. The novel fraction and the divergence series are also measured here |
| R3-major-1 | does this work on genomes other than human? | [`a3-cross-genome`](a3-cross-genome) | Arabidopsis Col-CEN and Ey15-2, label composition per feature set |
| R2/R3 smoothing | what does smoothing actually do, and to how much? | [`a4-smoothing-accounting`](a4-smoothing-accounting) | novel and multilabel accounting per feature set. 2.5643% of positions are labelled `novel` before smoothing, identical across all six sets |
| R3 divergence | how far from the reference does it still work? | [`a5-divergence-simulation`](a5-divergence-simulation) | 13 mutated-genome conditions, 100% coverage |
| Minor 1 | a parameter table with justifications, and a sensitivity note | [`a6-parameter-sweeps`](a6-parameter-sweeps) | five sweeps, and **no published number moves** |
| 2026-09-07 review | four factual checks and one mechanism check before resubmission | [`a7-review-checks`](a7-review-checks) | Ey15-2 is **18.66 % unrepresented** before smoothing (the 105 bp was the imputation residue). The gene-set percentages in the Results are **label-permuted** (exon is 11.5 %, not 49.4 %). Smoothing is **exactly** reverse-complement symmetric on both backends. The acrocentric flank rule holds once "flank" means the nearest label on the same branch. D4Z4 same-contig is by construction |
| — | the concordance harness used by the analyses above | [`concordance`](concordance) | reproduces the manuscript exactly on both HG002 haplotypes |

## The parameter sweeps, in one table

All five are in [`a6-parameter-sweeps`](a6-parameter-sweeps). Each has its own
`RESULTS.md` with the full numbers and a verdict for the supplementary table.

| sweep | parameter(s) | result |
|---|---|---|
| [`maxgap`](a6-parameter-sweeps/maxgap) | smoothing `max_gap` = 1000 bp | concordance is **flat to 0.0136 points** from 0 to 10,000 bp, so concordance could not have selected the value. It is justified instead from the distribution of gap lengths (99.4% of novel bases in runs <= 1 kb, median novel run exactly 31 bp = *k*) |
| [`tar1-its`](a6-parameter-sweeps/tar1-its) | TAR1/ITS merge distance = 1 kb | **no knee** in fragment count. The evidence is the bridging ceiling, which is stable up to d = 2000 and then doubles, so 1 kb sits in a wide safe range |
| [`rob`](a6-parameter-sweeps/rob) (result table withheld, see below) | Robertsonian (ROB) window / coverage / dominance | the three parameters named in the response letter control the ISCN *interval*, not the breakpoint: **the published breakpoint is identical across all 39 configurations** (3 samples × 13 settings: window 25–500 kb, minimum window fraction 0.1–0.5, specificity 0.5–0.9, each including the default). The two that do reach it move it by at most 9.1 kb over a 100-fold range |
| [`centromere-clustering`](a6-parameter-sweeps/centromere-clustering) | `--sil-threshold` 0.5, `--centroid-sd` 5 | **13 chromosomes retain a split at all 30 settings**. The four chromosomes in the allele heatmap are invariant at all 30, and 0.5 and 0.6 split the same 15 chromosomes |
| [`d4z4`](a6-parameter-sweeps/d4z4) | ten merge thresholds | **the published repeat count is 28,580 in all 51 runs**, identical under every threshold. Five of the ten thresholds have no effect, and the one sensitive parameter (the minimum array element size) is reported in Supplementary Table 3 |

## Conventions used throughout

Read these before re-running anything.

* **All runtimes were measured on the same hardware**, one Intel Xeon Gold 6342
  node (Ice Lake, 2.80 GHz), as in Supplementary Table 1. Jobs pin it with
  `--constraint=2.80GHz&Ice_Lake`.
* **Coordinates are stated in one frame.** Scaffolded and original-assembly
  coordinates differ on any flipped contig, and mixing them gives a plausible but
  wrong result, off by megabases, rather than a crash. Each table in the A2
  write-ups names its frame.
* **Never glob across the presmoothed/smoothed pair.** `*.smoothed.bed.gz` also
  matches `presmoothed`, and it sorts *first*, so `ls … | head -1` silently
  selects the wrong input. Scripts here name files exactly.
* **A control that must reproduce a known answer runs in every sweep.** The *k*
  sweep runs *k* = 31 through the same driver and requires it to reproduce the
  shipped pipeline to the base. The centromere sweep runs the published (0.5, 5)
  through all of steps 4–7. A sweep whose control fails is discarded, not
  interpreted.

## External code

Three analyses call code from other repositories, recorded below by repository
and commit. Copies of the first two are in `../centromere/scripts/` and
`../d4z4/`:

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

Each analysis directory contains `scripts/` (as run, including SLURM directives) and
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
  ROB sweep (`a6-parameter-sweeps/rob/`), and A7 Job 3 (`a7-review-checks/job3/`)
  read every ROB scaffold name, coordinate, and sample-specific value from
  `$ROBS_CONFIG`, sections `a2_step3`, `a6_rob`, and `a7_job3` of
  `../robs/conf/config.template.yaml`. At `b7d2092` these values were written in
  the scripts. With the filled configuration (dbGaP phs003920) each script gives
  output identical to the original's on the same inputs. `a2-k-sensitivity/step3/README.md`
  replaces the Step 3 part of `a2-k-sensitivity/README.md`, and
  `a6-parameter-sweeps/rob/README.md` summarises `rob/RESULTS.md` without
  coordinates.
- **Withheld pending the dbGaP ruling on the ROB assemblies:** the result
  tables and write-ups that contain ROB coordinates or scaffold names:
  `a2-k-sensitivity/results/sst1_breakpoints_by_k.tsv`, the Step 3 part of
  `a2-k-sensitivity/README.md`, `a6-parameter-sweeps/rob/results/rob_param_sweep.tsv`
  and `rob/RESULTS.md`, the A7 Job 3 provenance logs, and the ROB rows of
  `a7-review-checks/results/rc_symmetry_*.tsv` (the HG002 rows are complete).
  The Job 3 paragraph in `a7-review-checks/RESULTS.md` does not name the fusion
  scaffold.
- `a6-parameter-sweeps/d4z4/base_dir.patch` is dropped (see External code).
- Cluster paths are replaced by the variables in the repository README (Paths).
  Result tables that recorded input paths (`a2-k-sensitivity/step4/results/*_windows.tsv`)
  have the same substitution.
