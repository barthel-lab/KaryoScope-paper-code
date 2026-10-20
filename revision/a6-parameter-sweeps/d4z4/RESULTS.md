# A6 Sweep 4 — D4Z4 merge thresholds

> ## ⚠ CORRECTION (2026-09-03) — `BSR_MERGE_THRESHOLD` is NOT inert
>
> **§2 below wrongly lists `BSR_MERGE_THRESHOLD` among the "completely inert"
> constants and §4 wrongly recommends not disclosing it. Both are artefacts of how
> this sweep overrode the constant.**
>
> `analyze_d4z4_arrays.py` binds four constants as **default arguments**, evaluated
> at `def` time — `get_bsr_units` (:81, called bare at :1763), `get_bsr_clusters`
> (:990), `count_bsr_clusters` (:1664) all default to `BSR_MERGE_THRESHOLD`, and
> `detect_multiple_tracts` (:999) to `TRACT_GAP_THRESHOLD`. `sweep_d4z4.py` sets the
> module attribute, which **never reaches them**, so those two constants were never
> actually varied. The sweep still ran and still reported a flat line.
>
> Re-measured with `f.__defaults__` patched too (`scripts/sweep_corrected_default_binding.py`,
> results in `results/bsr_merge_and_tract_gap_CORRECTED.tsv`), on the same published inputs:
>
> | `BSR_MERGE_THRESHOLD` | 200 | 300 | 350–375 | **400–450** | 500–550 | 600 | 800 | 1200 |
> |---|---|---|---|---|---|---|---|---|
> | total repeat units | 37,724 | 30,521 | 30,521 | **28,580** | 21,275 | 19,843 | 16,993 | 3,632 |
> | HMM Pearson *r* | 0.893 | 0.936 | 0.936 | **0.9977** | 0.996 | 0.994 | 0.997 | 0.328 |
>
> **Total repeat units vary 2.2-fold across 200–800, the exact range called inert
> below.** The published 400 sits on a narrow 400–450 plateau and is a genuine local
> optimum for the external HMM correlation (0.9977 against 0.936 immediately below).
>
> **Consequences:** (a) `BSR_MERGE_THRESHOLD = 400` **must be disclosed** in the
> parameter table, with its real margin; (b) the "nine of the fourteen thresholds are
> completely inert" count in §2 drops to **eight**; (c) `TRACT_GAP_THRESHOLD` was
> re-tested *with* the fix and **is** genuinely inert (identical output at
> 5,000–20,000), so that row stands; (d) the other twelve constants are read as
> globals at call time and are unaffected — every other conclusion here holds.
>
> Found while testing whether the D4Z4 census could be salvaged at other *k*; see
> `../../a2-k-sensitivity/step4/RESULTS.md`.

Jobs 37146064 and 37146086, Ice Lake, **16 min 21 s** for **71 points over all 14 constants** (the second run extends the first to the four remaining undisclosed constants, so the inventory is fully covered). Floris's
`analyze_d4z4_arrays.py` (`barthel-lab/HPRC-D4Z4`) runs with its logic
**unmodified** — every swept constant is a module-level global read at call time,
so the driver overrides the module attribute and calls `main()`. The only source
edit is one line making `base_dir` env-overridable, because the upstream value is
a macOS path; carried as `base_dir.patch`.

Scope: chr4 and chr10, the only two chromosomes the script processes, **582
haplotypes**. Baseline at all published values: 629 arrays, 621 canonical,
8 degraded, 14,433 main-array repeat units, **28,580 total repeat units**.

Deliverable `d4z4_param_sweep.tsv` (71 rows). The `MIN_ARRAY_ELEMENT_SIZE`
refinement (2026-08-27, jobs 37157759 / 37157867) adds
`min_array_element_size_refined.tsv` (24 rows), `array_count_edges.tsv`,
`plateau.json` and `qualifying_sizes_full.json`, with drivers
`refine_min_array_element.py` and `confirm_array_edge.py`.

---

## 1. Headline — the published repeat count does not move at all

**`repeat_count_total` is 28,580 in all 71 runs.** Every one of the fourteen
thresholds, swept over its full range, leaves the total D4Z4 repeat-unit count
bit-identical. That is the published quantity, so this is the answer to the
sensitivity question in its strongest form.

The main-array repeat count (`repeat_count_main`, 14,433 at baseline) does move,
but only at range edges: 13,428–14,550 across all 71 runs, a −7.0%/+0.8% span,
and unchanged at the published value and its immediate neighbours for every
parameter.

## 2. Nine of the fourteen thresholds are completely inert

Zero haplotypes change, and array and repeat counts are bit-identical, over the
entire swept range:

| constant | published | swept | group |
|---|---|---|---|
| `NOVEL_MERGE_THRESHOLD` | 0.7 | 0.5 – 0.9 | 7a, stated |
| `EXON_MERGE_THRESHOLD` | 0.6 | 0.4 – 0.8 | 7a, stated |
| `ARRAY_GAP_THRESHOLD` | 10,000 | 5,000 – 20,000 | 7a, stated |
| `EXON_DENSE_THRESHOLD` | 0.5 | 0.3 – 0.7 | 7b, undisclosed |
| `BSR_MERGE_THRESHOLD` | 400 | 200 – 800 | 7c, undisclosed |
| `PLAM_MAX_EXON_FRACTION` | 0.35 | 0.25 – 0.45 | 7c, undisclosed |
| `WINDOW_SIZE` | 250,000 | 150,000 – 500,000 | 7c, undisclosed |
| `TRUNCATION_THRESHOLD` | 1,000 | 500 – 2,000 | 7c, undisclosed |
| `TRACT_GAP_THRESHOLD` | 10,000 | 5,000 – 20,000 | 7c, undisclosed |

**This is a measurement, not a failure to patch.** The control is that the other
five constants *do* move results through the same mechanism — `MIN_ARRAY_ELEMENT_SIZE`
at 2,000 changes 45 haplotypes and drops the array count from 629 to 572. So the
override works and the inert ones are genuinely inert on this dataset.

Mechanistically this is expected for the two gap-composition thresholds: the merge
path is entered only when two arrays sit within `ARRAY_MERGE_GAP`, which is rare
(552 of 582 haplotypes carry exactly one array), and the gaps that do occur are
either almost entirely novel or almost entirely not — so a threshold anywhere in
0.5–0.9 sits in a dead zone.

## 3. The five that do something, and how much margin the published value has

Each moves only at the edge of its range, and always on one side. (Ranges here
are at the swept grid's resolution except `MIN_ARRAY_ELEMENT_SIZE`, whose edges
were resolved exactly.)

<!-- MIN_ARRAY_ELEMENT_SIZE is split across two rows: its plateau differs by
     which count you read. Both edges resolved to the base pair, not to a grid
     step -- see the note below the table. -->

| constant | published | flat interval containing it | first setting that moves | effect there |
|---|---|---|---|---|
| `ARRAY_MERGE_GAP` | 4,000 | 2,000 – 6,000 | 1,000 / 8,000 | 1 hap / 4 haps |
| `COMBINED_MERGE_THRESHOLD` | 0.8 | 0.6 – 0.8 | 0.9 | 3 haps (22 at 1.0) |
| `BOUNDARY_WINDOW_SIZE` | 3,000 | 2,000 – 5,000 | 1,000 | 8 haps |
| `BOUNDARY_BLOCK_THRESHOLD` | 0.3 | 0.2 – 0.5 | 0.1 | 12 haps |
| **`MIN_ARRAY_ELEMENT_SIZE`** (array counts) | **1,500** | **901 – 1,600** | 1,601 | 5 haps, arrays 629 → 616 |
| **`MIN_ARRAY_ELEMENT_SIZE`** (`repeat_count_main`) | **1,500** | **1,201 – 1,500** | 1,501 | 2 repeat units, 0 haps |

Two of the brief's specific worries are now settled:

* **`ARRAY_MERGE_GAP` = 4 kb holds.** The manuscript defends it mechanistically as
  "a buffered cutoff above the canonical ~3.3 kb D4Z4 repeat-unit length", and it
  has ±2 kb of margin before anything changes. The mechanistic justification is
  supported, which is worth more than a fitted optimum.
* **`BOUNDARY_BLOCK_THRESHOLD` = 0.3 is not exposed after all.** The brief called
  it the most exposed because the manuscript describes it as purely empirical.
  Measured, it has 0.1 of margin below and 0.2 above. It is fine.

`MIN_ARRAY_ELEMENT_SIZE` is the outlier, and its margin depends on which count
you read. Array counts and the canonical/degraded split are invariant over
**901–1,600**, so the published 1,500 is interior with 100 bp above it. The
main-array repeat count is invariant over exactly **1,201–1,500**: the published
value is the **upper edge**, because 1,500 bp is itself a populated cluster size
and the gate is `size >=`. Crossing it costs 2 repeat units of 14,433 and no
haplotype changes its array count. It is the most sensitive of the fourteen by a
wide margin.

Both edges were resolved to the base pair rather than to a grid step. The
step-250 grid could only say the first change lay somewhere in (1500, 1750]; an
earlier version of this document reported that 250 bp as if it were the margin,
which it was not. Because the presmoothed input is binned at 100 bp, **all 33
qualifying cluster sizes are multiples of 100**, so 100 bp is the finest
meaningful step for this parameter and no finer grid can move the boundary. Full
method and per-setting table:
`docs/write-ups/a6-d4z4-min-array-element-plateau-REFINEMENT.md` (jobs 37157759,
37157867).

## 4. Disclosure recommendation

| constant | recommendation |
|---|---|
| **`MIN_ARRAY_ELEMENT_SIZE` = 1,500** | **Disclose.** It is the most sensitive parameter of the fourteen. The published value is the exact upper edge of the interval over which the main-array repeat count is invariant — 1,501 moves it by 2 units of 14,433 — though array counts hold to 1,600 and the published *total* repeat count is inert across the whole swept range. There is a second reason: it coincidentally shares its value with the manuscript's 1,500 bp *boundary-extension distance*, which is a different quantity — disclosing it with its real meaning ("minimum size for a full array element, used for counting") removes a genuine ambiguity rather than adding a row. |
| `BSR_MERGE_THRESHOLD` = 400 | Do not disclose. Inert over 200–800; implementation detail. Its in-code comment calls it "conservative", which the measurement supports. (The `get_bsr_units()` docstring says "default 500bp"; the code is 400 and 400 is what ran. The mismatch is real but, given the parameter is inert over 200–800, it changes nothing.) |
| `PLAM_MAX_EXON_FRACTION` = 0.35, `WINDOW_SIZE` = 250,000, `TRUNCATION_THRESHOLD` = 1,000, `TRACT_GAP_THRESHOLD` = 10,000 | Do not disclose. All four inert over the ranges above. Swept only so the table can say what was measured rather than what was not looked at. |
| `EXON_DENSE_THRESHOLD` = 0.5 | Do not disclose as a tuned parameter, but note it. It defines the primitive the array definition rests on, so it *sounds* load-bearing; it is inert over 0.3–0.7. Worth one sentence saying so, since a reader who finds it will otherwise wonder. |
| the seven stated thresholds | Keep all seven, but **write the justifications to say what was measured.** Three of them (`NOVEL_MERGE_THRESHOLD`, `EXON_MERGE_THRESHOLD`, `ARRAY_GAP_THRESHOLD`) have no measurable effect on this dataset, and a justification that implies they were tuned would be wrong. "Stated for reproducibility; no effect on either published count across [range]" is both honest and stronger. |

## 5. Caveats

* **582 haplotypes, not 231.** The manuscript's *n* = 231 is the subset with
  external HMM-based repeat counts (the correlation in `plot_hmm_correlation.py`,
  a separate reviewer comment). This sweep covers the full chr4 + chr10 set the
  array caller processes.
* **`haps_changed` counts haplotypes whose *array count* changed**, so it can read
  0 while repeat counts shift slightly — `MIN_ARRAY_ELEMENT_SIZE` = 1,000 moves
  `repeat_count_main` by 4 with no array-count change. Read the count columns too.
* **One at a time only.** Interactions were not swept. Given that nine constants
  are inert and the rest move only at range edges, a full factorial is unlikely to
  be informative, but it was not measured.

Related: [[a6-maxgap-sweep-result]] [[a6-inputs-found-2026-08-26]]
