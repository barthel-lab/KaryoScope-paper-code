# A2 Step 3: the Robertsonian breakpoint as a function of *k*

Step 3 re-derives the published Robertsonian (ROB) breakpoint from chromosome
annotations made with each of the six *k*-length databases of Step 2
(`HKS_human_CHM13_v2_k{11,21,31,41,51,61}`), on the three ROB cell-line
assemblies (GM03417, GM03786, GM04890).

## What it backs in the manuscript

- **Methods**, Robertsonian breakpoint refinement: "We repeated the Robertsonian
  breakpoint refinement using chromosome annotations from each of the six
  databases ... Across *k* = 21-61, the inferred breakpoint moved by 900 bp
  (GM03417), 50 bp (GM03786), and 0 bp (GM04890). At *k* = 11, refinement
  failed on two of the three samples ...".
- **Results**, ROB section: "moved the inferred breakpoint by at most 900 bp
  (see Methods)".
- **Reporting Summary**, replication: "Robertsonian breakpoints moved by at
  most 900 bp".

The per-*k* result table (`sst1_breakpoints_by_k.tsv`) holds ROB coordinates
and is not distributed; a run of these scripts on the dbGaP assemblies
regenerates it.

**Provenance of the table.** The *k* = 31 rows were re-run on 2026-10-20 and are
identical, field for field, to the published table for all three samples. The
*k* = 11, 21, 41, 51 and 61 rows are from the original 2026-08 runs; those
databases were not retained. They were built on 2026-08-03 with `karyoscope
build` at KaryoScope `9ccba0b` and hks `c5fecf8` (see the repository README,
Versions), from the six CHM13 v2 feature-set BEDs and priority files with a
fixed *k* per database. The spec generator used then was not kept. The public
recipe for the shipped database, KaryoScope `docs/recipes/human-chm13v2.md`
(`kmer: s: 31`), builds the same feature sets; setting `s` to 11, 21, 41, 51
or 61 gives the other databases. A database built that way has not been
compared with the originals.

## Reproduction path

1. `run_step1_scaffold.sh` (public KaryoScope): `karyoscope annotate` (region and
   chromosome sets) and `karyoscope scaffold --input hap1=<fasta> --feature-set
   region --mode bed` on haplotype 1 of each sample. This writes the *k* = 31
   scaffold map that step 3 holds fixed.
2. `run_step3_sweep.sh` + `sweep_refine.py` per (sample, *k*), with
   `$KS_ISCN_SRC` pointing at `../../../vendored/KaryoScope-ISCN/src` (the
   KaryoScope-ISCN analysis modules at `5507369`).

The original runs made the scaffold map with `karyoscope-iscn analyze`
(`run_step1_baseline.sh`, below). `analyze` builds it by calling the same
`karyoscope scaffold` command, and the map is the only `analyze` output step 3
uses: on 2026-10-20, `run_step1_scaffold.sh` (KaryoScope 2.3.0) gave a scaffold
map byte-identical to `analyze`'s on all three samples, both at KaryoScope-ISCN
`5507369` and at its current `main`. With that map, the vendored modules and the
ROB configuration, the *k* = 31 rows above are reproduced.

## Inputs and configuration

- The ROB assemblies (dbGaP phs003920; see `../../../robs/README.md`, step 1),
  under `$KS_GUARRACINO/ks-robs/results/<sample>/assembly/haplotype{1,2}/`.
- The shipped databases `HKS_human_CHM13_v2`, `HKS_human_CHM13_cytoband` and
  `KS_human_CHM13` (v1) in the KaryoScope database registry, and the six
  *k*-length databases of Step 2 under
  `$KS_WORK/ks_build_sources/chm13v2_ksweep/db` (see the repository README,
  Versions).
- The KaryoScope-ISCN analysis modules vendored in this repository
  (`../../../vendored/KaryoScope-ISCN/src`, commit `5507369`), via
  `$KS_ISCN_SRC`. The as-run `analyze` scripts (steps 1, 1b, 2) also need the
  `karyoscope-iscn` package itself, which is not public; they are not on the
  reproduction path.
- **The ROB configuration** (`$ROBS_CONFIG`): the `a2_step3` section of
  `../../../robs/conf/config.template.yaml`, filled from your own annotation.
  Per sample it names the `*_robertsonian` breakpoint group (fusion scaffold,
  ISCN window, partner chromosomes) and gives the *k* = 31 detected position,
  the scaffold length, the SST1 block midpoint and the comparison values. No
  coordinate or scaffold name is written in the scripts.

Set `KS_WORK`, `KS_SOFTWARE`, `KS_GUARRACINO`, `KS_ISCN_SRC` and `ROBS_CONFIG`, and copy
`scripts/*.py` into `$KS_WORK/a2_step3` (the drivers call `sweep_refine.py` and
`binsweep_analyse.py` there, and the Python scripts read the configuration
through `rob_samples.py` next to them).

## Scripts, in run order

| step | script | does |
|---|---|---|
| 1a | `run_step1_scaffold.sh` (SLURM array 1-3) | **reproduction path**: the *k* = 31 scaffold map with public `karyoscope annotate` + `karyoscope scaffold` (see above) |
| 1 | `run_step1_baseline.sh` (SLURM array 1-3) | as run: `karyoscope-iscn analyze` (KaryoScope-ISCN `5507369`, not public) at *k* = 31 on both haplotypes; finds the fusion scaffold and the detected breakpoint; leaves the cached intermediates in `$KS_WORK/a2_step3/<sample>/work` |
| 1b | `run_step1_verify.sh` | re-runs the report on the cached intermediates after two ISCN false-positive fixes |
| 2 | `run_step3_gate.sh` | the same with `--refine-pivot chromosome`; checks the inferred breakpoint against the SST1 midpoint |
| 2b | `sst1_blocks.sh` | SST1 blocks on each fusion scaffold, from the acrocentric annotation (gives `sst1_mid`) |
| 3 | `run_step3_sweep.sh` (array 1-18) + `sweep_refine.py` | per (sample, *k*): annotate the chromosome set, clip to the *k* = 31 scaffold map, remap onto that fixed layout, run the shipped `refine_breakpoint` + `cumulative_breakpoint`; one row per task (`row.tsv`), concatenated into `sst1_breakpoints_by_k.tsv`. *k* = 31 is the control: it must reproduce step 2 to the base |
| 4 | `run_step3_sweep_binned.sh` (array 1-18), `run_binsweep_analyse.sh` + `binsweep_analyse.py` | the same sweep on 100 bp winner-take-all binned input, the representation of the submitted manuscript |
| 5 | `run_presmoothed_ci_test.sh`, `run_binned100_test.sh`, `binned100_analyse.py`, `ci_containment.py`, `reproduce_published_density.py`, `score.sh`, `gap.sh` | reconciliation of the submitted 200-525 kb range and 99% interval widths: presmoothed vs smoothed input, binned vs per-base input, an independent implementation of the cumulative method over the whole scaffold, and the width of the unassignable gap |
| 6 | `run_v1_delta_test.sh`, `v1_vs_v2_compare.py` | the chromosome annotation with the v1 database `KS_human_CHM13`, compared with v2 |

`rob_samples.py` reads the configuration (`rob_samples.py sample scaffold sst1_mid`
prints one line per sample).

## Result (relative values only)

Across *k* = 21-61 the inferred breakpoint moves by at most 900 bp (900 / 50 /
0 bp for the three samples), against distances of 46-252 kb to the SST1
midpoint. At *k* = 11 refinement falls back on two of the three samples. The
submitted "approximately 200-525 kb" came from 100 bp binned input; on per-base
input the distances are 46-252 kb, and in the binned configuration the
breakpoint is not stable across *k*, so the manuscript uses the per-base
numbers.

## What is held fixed

Only the chromosome-identity evidence fed to refinement varies. Held fixed: the
*k* = 31 scaffold layout (a per-*k* layout could flip a contig), the detected
event, and the database version.

Four traps met on the way, each of which gave a plausible wrong answer:

1. `ls *chromosome*smoothed*.bed.gz | head -1` selects `presmoothed` (it sorts
   first). The drivers name every file exactly.
2. A BED's extent depends on *k*: *k*-mer starts run `0..L-k`, so a BED at *k*
   reaches `31-k` bp past the *k* = 31 map's length and `remap-bed` refuses it.
   The drivers clip to the map's lengths (at most 30 bp at contig ends).
3. `iscn_regions.tsv` is in original-assembly coordinates, while the
   configuration and the scaffolded BEDs are scaffolded. One fusion contig is
   reverse-complemented in the scaffold, so mixing the frames gives a spurious
   error of several megabases.
4. Both haplotypes must be passed per sample: one haplotype is reported against
   a haploid baseline, which misreads these XX samples.

## Changes from KaryoScope-revisions `b7d2092`

- Every per-sample coordinate, scaffold name and comparison value is read from
  `$ROBS_CONFIG` (`a2_step3`); `rob_samples.py` is new.
- Cluster paths are replaced by the variables in the repository README (Paths);
  `#SBATCH --output` writes to the submission directory.
- Comments that quoted per-sample coordinates or values now point at the
  configuration keys.

Checked 2026-10-20 with the real configuration, at *k* = 31 (the other *k*
databases were not retained): the original scripts and the refactored scripts
were run in sequence on the same inputs, and every output file of every step
compares as follows.

| script | result |
|---|---|
| `run_step1_baseline.sh`, `run_step1_verify.sh`, `run_step3_gate.sh` | identical, except `iscn_segments.tsv`, whose `arms` column prints a Python `frozenset`; its element order varies between any two runs (string hashing), and the files are identical once that is normalised. No result uses the column |
| `run_presmoothed_ci_test.sh`, `run_binned100_test.sh`, `run_v1_delta_test.sh` | identical |
| `run_step3_sweep.sh` + `sweep_refine.py` | identical; the *k* = 31 rows equal the published table |
| `run_step3_sweep_binned.sh`, `run_binsweep_analyse.sh` + `binsweep_analyse.py` | identical |
| `score.sh`, `gap.sh`, `sst1_blocks.sh`, `ci_containment.py`, `reproduce_published_density.py`, `v1_vs_v2_compare.py`, `binned100_analyse.py` | identical (stdout and stderr) |
| `run_step1_scaffold.sh` (new) | scaffold maps identical to `karyoscope-iscn analyze`'s, at ISCN `5507369` and current `main` |

`sweep_refine.py` was also run against the vendored KaryoScope-ISCN modules
(`5507369`), and with the scaffold map from `analyze` at `5507369`: the *k* = 31 rows
are identical to the published table. The comparison above ran with the
KaryoScope-ISCN installed on the cluster at the time (`main`, `4eb7128`) for both
the original and refactored scripts.
