# A2: sensitivity to *k*

This directory holds Steps 0-4 of the *k* sweep.

- `scripts/run_step01.sh`: Step 0 (pre-flight) and Step 1, the k-mer uniqueness
  sweep on the `repeat` feature set (`hks node-stats`).
- `scripts/run_a2_step2.sh`, `run_a2_step2_k11.sh`: Step 2, annotation of each
  assembly at each *k* on the HKS backend, with per-base leaf, internal-node and
  novel fractions per feature set.
- `scripts/run_a2_concord.sh`: Step 2b, per-base concordance of each *k* against
  *k* = 31 in both directions.
- `step4/`: novel fraction, divergence, the D4Z4 census and the centromere
  clustering as a function of *k*. Read `step4/RESULTS.md`.

- `step3/`: the Robertsonian breakpoint re-derived at every *k*. The scripts
  read every per-sample scaffold name and coordinate from the ROB configuration
  (`$ROBS_CONFIG`, see `../../robs/README.md`), and run with dbGaP access to the
  ROB assemblies. Its result table is withheld (it holds ROB coordinates); the
  result, that the breakpoint moves at most 900 bp across *k* = 21-61, is stated
  in the manuscript. Read `step3/README.md`.
