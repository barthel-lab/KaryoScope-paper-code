# A6 ROB sweep: the Robertsonian breakpoint refinement parameters

`scripts/sweep_rob.py` (driver `scripts/run_rob_sweep.sh`, about a minute)
varies the parameters of the Robertsonian breakpoint refinement one at a time
and records, per sample and setting, the content-dominance interval
(`refine_breakpoint`) and the inferred breakpoint (`cumulative_breakpoint`).
It reads the *k* = 31 scaffolded chromosome BEDs written by A2 Step 3
(`$KS_WORK/a2_step3/sweep/<sample>/k31/chromosome.scaffolded.bed.gz`, from
`../../a2-k-sensitivity/step3/scripts/run_step3_sweep.sh`) and does not
re-annotate.

Three sweeps:

- **A**, `refine_breakpoint`: `window` (25-500 kb, default 100 kb),
  `min_window_fraction` (0.10-0.50, default 0.25), and `specificity` (0.50-0.90,
  default 0.50), with `cumulative_breakpoint` over the whole scaffold as shipped.
- **B**, `cumulative_breakpoint`: `resolution` (100 bp-10 kb, default 100 bp) and
  `plateau_threshold` (0.90-1.0, default 0.99).
- **C**, the A settings with `cumulative_breakpoint` over the refine interval
  +/- 1 Mb instead of the whole scaffold.

Held fixed: *k* = 31, the *k* = 31 scaffold layout, the detected event, the
database version, and the chromosome pivot.

## What it backs in the manuscript

All five rows of the "Robertsonian breakpoint refinement" block of
Supplementary Table 2 (`tab:S3_parameters_core`):

| row | value | evidence quoted | sweep |
|---|---|---|---|
| chromosome-assignment window | 100 kb | breakpoint unchanged over 25-500 kb. Transition-interval centre moves up to 660 kb. At 25 kb refinement falls back to an unrefined interval | A, `window` |
| coverage threshold | 25% | unchanged over 10-50% | A, `min_window_fraction` |
| dominance threshold | 50% | unchanged over 50-90% | A, `specificity` |
| cumulative-curve grid | 100 bp | coarsening to 10 kb moves the breakpoint by at most 9.1 kb, far below the 46-252 kb distances to the SST1 midpoint | B, `resolution` |
| plateau threshold | 99% | identical over 90-100%, and plateau width moves from 30.2 Mb to 3.3 kb | B, `plateau_threshold` |

The grid row comes from this sweep (`resolution` in sweep B), not from the
A2 Step 3 binned analyses. The result table (`rob_param_sweep.tsv`) has 108 rows from three sweeps. Sweep A (39 rows) varies the window,
minimum window fraction, and specificity, 13 settings per sample. Sweep B (30 rows)
varies the cumulative-curve grid and the plateau threshold, 10 settings per sample.
Sweep C (39 rows) repeats sweep A's settings with the earlier 1 Mb cumulative window. The table
contains Robertsonian (ROB) coordinates and is not distributed.

## Run

```bash
export KS_WORK=... KS_SOFTWARE=... ROBS_CONFIG=<filled ROB configuration>
export KS_ISCN_SRC=<this repository>/vendored/KaryoScope-ISCN/src   # ISCN modules at 5507369
mkdir -p $KS_WORK/a6_sweeps/rob && cp scripts/sweep_rob.py $KS_WORK/a6_sweeps/rob/
sbatch scripts/run_rob_sweep.sh      # writes $KS_WORK/a6_sweeps/rob/rob_param_sweep.tsv
```

The per-sample values (fusion scaffold, partner chromosomes, *k* = 31 detected
position, scaffold length, SST1 midpoint) come from the `a2_step3` section of
the configuration, named in the `a6_rob` section
(`../../../robs/conf/config.template.yaml`).

## Changes from KaryoScope-revisions `b7d2092`

The per-sample values are read from `$ROBS_CONFIG` instead of being written in
the script. Cluster paths are replaced by the variables in the repository
README (Paths). KaryoScope-ISCN is imported from the vendored copy at `5507369`
(`$KS_ISCN_SRC`). The sweep fails with KaryoScope-ISCN versions after 2026-09-11,
which dropped `DEFAULT_WINDOW_PAD`.

Checked 2026-10-20 with the real configuration: the original script (with
KaryoScope-ISCN `5507369`) and the refactored script (with the vendored modules)
both write a `rob_param_sweep.tsv` identical to the published table (108 rows),
on sweep inputs regenerated that day. Running `run_rob_sweep.sh` end to end gives
the same result.
