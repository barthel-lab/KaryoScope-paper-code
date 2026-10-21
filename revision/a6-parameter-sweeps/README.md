# A6 — parameter sweeps

Five sweeps behind the supplementary parameter table. Each subdirectory has its
own `RESULTS.md` with the numbers and a verdict.

The response letter promised "the compact parameter table with justifications,
and a brief sensitivity note for the smoothing gap, ROB window/dominance, D4Z4
thresholds, and centromere clustering criteria". Every row needs a value, where
it is used, why that value, and **what evidence backs it** — that last column is
what these sweeps supply.

| sweep | cost | verdict |
|---|---|---|
| [`maxgap`](maxgap) | 9–28 s per point | value fine, **stated reason was wrong** — rewritten to justify from the gap-length distribution |
| [`tar1-its`](tar1-its) | minutes | no knee; the bridging ceiling is the evidence |
| [`rob`](rob) (scripts; result table withheld) | **44 s total** | the brief named the wrong parameters; the published breakpoint is invariant |
| [`centromere-clustering`](centromere-clustering) | **84 s total** | structure survives; one open question on chr19 |
| [`d4z4`](d4z4) | ~16 s per point | published repeat count invariant across all 51 runs; 5 of 10 thresholds inert; recommends disclosing `MIN_ARRAY_ELEMENT_SIZE` |

## Two standing rules these follow

**Do not tune a threshold to make the sweep look good, and do not quietly drop a
bad row.** A sweep that shows sensitivity is a legitimate outcome to report. Of
the five, one did change a manuscript sentence (`maxgap` — the stated rationale
was unsupported) and none changed a published number.

**Sweeps are cheap because they run on existing intermediates.** `maxgap` smooths
an existing presmoothed BED; `rob` reads the existing *k* = 31 scaffolded BEDs;
`centromere-clustering` re-gates one clustering run. None of them re-annotates.
If a sweep design implies re-annotation, check first whether the parameter is
consumed downstream of the expensive step — three of these five were.
