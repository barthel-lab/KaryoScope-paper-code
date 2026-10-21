# A6 — parameter sweeps

These five sweeps support the supplementary parameter table. Each subdirectory
has its own `RESULTS.md` with the numbers and a verdict.

In the response letter we promised "the compact parameter table with justifications,
and a brief sensitivity note for the smoothing gap, ROB window/dominance, D4Z4
thresholds, and centromere clustering criteria". Each row of the table gives a
value, where it is used, why that value was chosen, and **the evidence for it**.
These sweeps provide the evidence.

| sweep | cost | verdict |
|---|---|---|
| [`maxgap`](maxgap) | 9–28 s per point | value fine, **stated reason was wrong**, rewritten to justify it from the distribution of gap lengths |
| [`tar1-its`](tar1-its) | minutes | no knee (the bridging ceiling is the evidence) |
| [`rob`](rob) (scripts only, result table withheld) | **44 s total** | the brief named the wrong parameters, and the published breakpoint is invariant |
| [`centromere-clustering`](centromere-clustering) | **84 s total** | structure survives. The FISH-validated chr19 variant stays in the major cluster because array lengths are not clustering features |
| [`d4z4`](d4z4) | ~16 s per point | published repeat count invariant across all 51 runs, 5 of 10 thresholds inert, and `MIN_ARRAY_ELEMENT_SIZE` is reported in Supplementary Table 3 |

## Two rules for the sweeps

**Do not tune a threshold to make the sweep look good, and do not quietly drop a
bad row.** A sweep that shows sensitivity is a result to report. One of the
five changed a manuscript sentence (`maxgap`, whose stated rationale was
unsupported), and none changed a published number.

**Sweeps are cheap because they run on existing intermediates.** `maxgap` smooths
an existing presmoothed BED, `rob` reads the existing *k* = 31 scaffolded BEDs,
and `centromere-clustering` re-gates one clustering run. None of them re-annotates.
If a sweep design seems to need re-annotation, first check whether the parameter
is used after the expensive step. For three of these five it was.
