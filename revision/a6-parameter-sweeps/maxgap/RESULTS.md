# A6 §3 — smoothing `max_gap`: results

**Run:** 2026-08-24, Intel Xeon Gold 6342 (Ice Lake) @ 2.80 GHz, the Table S1 node class.
**Jobs:** `36807085` (multilabel run lengths, array 1–12), `36807757` (max_gap concordance, array 1–8).
**Verdict:** **sweep-backed row — but the manuscript's stated rationale is wrong and should change.**

---

## 1. Headline

Two things were measured. They point the same way, and neither is what `main.tex:366`
currently claims.

1. **Concordance is flat in `max_gap`.** Across 0 → 10,000 bp the post-smoothing
   equivalent fraction moves from 95.2344% to 95.2480% — a total spread of
   **0.0136 percentage points** over four orders of magnitude.
2. **The gap-length distribution does discriminate**, and it is what actually
   justifies 1000 bp.

So the value is safe, but **the reason given for it is not the reason that supports it**.

## 2. Concordance sweep

HG002 hap1, `repeat` feature set, RepeatMasker ground truth. The `repeat` set was
used because it is the one whose truth is in hand *and* whose anchor is proven exact
(job 35385021 reproduced Supp Note 1 to the base). Annotation was not re-run:
`hks smooth` takes the existing presmoothed BED and `--max-gap`, so each point cost
9–28 s of smoothing plus one concordance pass.

Pre-smoothing baseline, identical in every run because it is the same input:

| class | bases | fraction |
|---|---|---|
| equivalent | 2,492,713,967 | 84.7787% |
| compatible | 271,820,734 | 9.2448% |
| novel | 75,396,753 | 2.5643% |
| incompatible | 100,328,647 | 3.4122% |

Post-smoothing, by `max_gap`:

| max_gap | equivalent | compatible | novel | incompatible |
|---|---|---|---|---|
| 0 | 95.2344% | 1.0757% | 0.0406% | 3.6493% |
| 100 | 95.2398% | 1.0700% | 0.0406% | 3.6496% |
| 250 | 95.2421% | 1.0677% | 0.0406% | 3.6496% |
| 500 | 95.2440% | 1.0657% | 0.0406% | 3.6496% |
| **1000** | **95.2461%** | **1.0637%** | **0.0406%** | **3.6496%** |
| 2000 | 95.2473% | 1.0624% | 0.0406% | 3.6496% |
| 5000 | 95.2479% | 1.0618% | 0.0406% | 3.6496% |
| 10000 | 95.2480% | 1.0618% | 0.0406% | 3.6496% |

`max_gap = 1000` reproduces the published 95.25% exactly.

The curve is monotone and saturates — 5000 and 10000 differ by 0.0001 points — so
there is a knee, and 1000 bp is past the steep part of it. But the whole dynamic
range is 0.014 points. **Concordance cannot distinguish 100 bp from 10 kb**, so it
cannot be the evidence that selected 1000.

### The control does not exist as the brief assumed

§3 item 3 expected `max_gap = 0` to be "smoothing switched off" and to reproduce the
pre-smoothing numbers. **It does not, and it should not.** `--max-gap` is the maximum
*coordinate gap between adjacent intervals considered connected*; at 0, touching
intervals are still connected, so hierarchy-aware smoothing still runs in full. That
is why `max_gap = 0` gives 95.2344%, not 84.7787%.

The smoothing-off baseline is the **presmoothed row**, which is already in every
summary. Read against it, the split is:

* smoothing itself: 84.7787% → 95.2344%, **+10.46 points**
* the entire `max_gap` range on top of that: **+0.014 points**

This is not a harness failure — it is the parameter meaning something narrower than
the brief assumed. Worth saying plainly in the table so nobody re-derives it.

## 3. Multilabel run lengths — the measurement that was missing

`main.tex:366` says novel k-mers "are treated uniformly with multilabel k-mers during
smoothing", so `max_gap` bridges gaps from both populations, but only the novel one
had been characterised (A4). Multilabel is defined here as *the presmoothed label is
an internal node of the feature set's hierarchy* — a k-mer matching several leaves
resolves to their LCA, so an internal-node call is the ambiguous call. Method mirrors
`run_a4_novel_v3.sh` exactly (same join, same N masking, same buckets), so the tables
are directly comparable.

### Multilabel is the larger population, by a lot

| set | multilabel pre | multilabel post | median resolved run |
|---|---|---|---|
| chromosome | **16.03%** | 1.72% | 8 bp |
| repeat | 9.59% | 0.99% | 8 bp |
| region | 8.74% | 0.85% | 5 bp |
| acrocentric | 0.08% | 0.04% | 3 bp |
| subtelomeric | 0.013% | 0.007% | 3 bp |
| gene | 0% | 0% | — |

(HG002 hap1; CHM13 is within a point of these except acrocentric.) Novel, for
comparison, is ~2.5% of the repeat set. **In the chromosome set 16% of evaluated
bases are ambiguous pre-smoothing, against 2.5% novel** — so the population the
manuscript mentions in passing is the dominant one, and smoothing takes it to 1.7%.

`gene` correctly reports zero: it is a flat three-leaf tiling with no internal nodes,
which is a useful check that the leaf/internal split is being read right.

### Fraction below 1 kb — the number that justifies the threshold

| set | HG002 hap1 | CHM13 |
|---|---|---|
| subtelomeric | 100.00% | 100.00% |
| region | 99.50% | 99.31% |
| repeat | 98.81% | 98.59% |
| chromosome | 95.56% | 93.12% |
| acrocentric | 93.25% | **67.31%** |
| *novel (A4, reference)* | *99.4%* | — |

Median resolved run is **3–9 bp** in every set — far below k=31, and three orders of
magnitude below the threshold.

### The acrocentric split is an N-masking artefact, not a CHM13 quirk

**CORRECTED 2026-08-24.** The earlier reading here — that 67.3% was "a
CHM13-specific feature" — was wrong.

CHM13 acrocentric sits at 67.3% below 1 kb against HG002 hap1's 93.25%. That gap
is not biology and not a CHM13 peculiarity: **HG002 cannot contribute the relevant
sequence at all.**

HG002 hap1's N runs are exactly the acrocentric short arms:

```
chr14_hap1  3,463,694-4,966,831
chr15_hap1  2,973,704-5,724,253
chr21_hap1  2,675,312-4,976,002
chr22_hap1  2,633,586-3,713,463      total 7.63 Mb
```

Following A4's method, N runs are masked out of the measurement — they are absence,
not ambiguity. That masking removes precisely the rDNA and p-arm regions that
generate long multilabel runs, so HG002's 93.25% is measured over an assembly with
those regions deleted.

CHM13 is T2T and gapless, so nothing is masked. Every one of its 427 long (>=1 kb)
acrocentric multilabel runs falls on an acrocentric chromosome:

| chromosome | runs | bases |
|---|---|---|
| chr13 | 154 | 412,172 |
| chr21 | 113 | 297,502 |
| chr15 | 102 | 276,240 |
| chr22 | 24 | 123,874 |
| chr14 | 34 | 74,221 |
| **total** | **427** | **1,184,009** |

That is the entire 1-10 kb bucket, exactly.

**So 67.3% is the honest number and 93.25% is the flattered one.** In acrocentric
rDNA and p-arm sequence — which only a T2T assembly can measure — multilabel runs
are much longer, and `max_gap = 1000` bridges only about two thirds of that span.
That is a real limitation of the threshold in exactly the regions this paper is
about, and it should be stated as one rather than as a quirk of one genome.

It does not change the recommended value. The concordance sweep is flat, and raising
`max_gap` far enough to bridge 10 kb runs would start joining genuinely distinct
arrays. But the caveat column should say the threshold is calibrated to the bulk of
the distribution and leaves long acrocentric-array ambiguity unresolved by design.

## 4. Recommended row, and the manuscript change

**Row verdict: sweep-backed.** 1000 bp is defensible and the sweep supports it.

**But `main.tex:366` should be rewritten.** It currently says the value "was chosen by
evaluating the concordance between the smoothed chromosome feature annotations of the
HG002 assembly and the corresponding reference-derived ground truth". Those numbers
now exist, and they show concordance is flat to 0.014 points from 100 bp to 10 kb.
Publishing that sentence alongside this sweep invites the obvious question.

Suggested justification, which the evidence does support:

> `max_gap` is set to 1000 bp because it sits just above the length distribution of
> the gaps it has to bridge. Novel k-mer runs are 99.4% shorter than 1 kb and
> multilabel runs 93–100% shorter, with a median resolved run of 3–9 bp. Concordance
> against the reference-derived ground truth is insensitive to the exact value —
> post-smoothing agreement varies by 0.014 percentage points between 100 bp and
> 10 kb — so the threshold is chosen from the gap distribution rather than tuned.

That is a stronger claim than the current one, because insensitivity is a better
answer to R3-major-3 than a fitted optimum would be.

## 5. Files

```
a6_sweeps/run_multilabel_runlen.sh, run_maxgap_sweep.sh
a6_sweeps/results/multilabel/multilabel_{asm}_{set}.tsv
a6_sweeps/results/multilabel/multilabel_runlen_hist_{asm}_{set}.tsv
a6_sweeps/results/maxgap/maxgap{G}_repeat.summary.tsv    G in 0..10000
```

## 6. Not done

* **Chromosome-set concordance.** The sweep ran on `repeat` because that truth is
  staged and anchored. `main.tex:366` names the *chromosome* set specifically, so
  matching it exactly needs a reference-derived chromosome truth built first.
  The multilabel measurement above *does* cover chromosome, and it is the part that
  carries the argument.
* TAR1/ITS (§6), centromere clustering (§5), D4Z4 (§7) — not started.
* ROB (§4) — **§4's premise is stale.** The 200–525 kb discrepancy was reconciled on
  2026-08-24 (binned input; addendum 6 of the delta brief) and option B (unbinned)
  was chosen. Supp. Fig. 5 has since been regenerated. ROB is blocked only on PR #4.
