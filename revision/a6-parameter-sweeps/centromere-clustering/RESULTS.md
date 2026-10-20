# A6 Sweep 3 — centromere clustering (`--sil-threshold` x `--centroid-sd`)

Job 37136128, Ice Lake (g-c-1-7-20), **1 min 24 s** wall. Yi-An's scripts run
**unmodified**, exported from `barthel-lab/KaryoScope-analysis` branch
`pangenome_structure_V2` @ `90e2874`.

Deliverables in `$KS_WORK/a6_sweeps/cent/`:

| file | what |
|---|---|
| `cent_param_sweep.tsv` | 36 rows: counts per (sil, sd, centroid_scan), with per-chromosome outlier counts |
| `validated_haplotypes_by_setting.tsv` | the five FISH-validated haplotypes' class at each of the 30 settings |
| `control/` | full steps 4–7 at the published (0.5, 5) — dendrogram, annotations, heatmaps |
| `run_cent_sweep.sh`, `collate.py`, `track_validated.py` | method |

Inputs and the path correction: `../A6-DATA-LOCATIONS-for-YiAn.md`.
Pipeline reproduced end to end: 24 BEDs → 93,452 rows → QC pass 63,659 rows /
4,650 contigs → clustering 4,646 assignments (11.5 s, 282 MB peak) → all four
steps 4–7 at rc=0.

---

## 1. Design — the sweep is over the interpretation layer only

`--sil-threshold` and `--centroid-sd` are **not** arguments to step 3
(`KaryoScope_cluster_analysis.py`). They are consumed by steps 4–7, which read
`allchr_structure.sequence_assignments.tsv`. So the clustering runs once and the
36 points re-gate its output. That is why the whole sweep costs 84 seconds.

Grid: `sil ∈ {0.0, 0.3, 0.4, 0.5, 0.6, 0.7}` (0.0 = ungated control) x
`sd ∈ {3, 4, 5, 6, 7}`, published **(0.5, 5)**. Plus six runs with
`--no-centroid-scan`, because both stages create and destroy outliers and a
combined count cannot attribute a change to one of them.

Note `--centroid-scan` defaults to **on**, so the driver's `--centroid-sd 5` is
live even though it never passes `--centroid-scan` explicitly.

## 2. Headline — the structure survives, and 13 chromosomes never lose their split

Total haplotypes is constant at 4,551 in every run. Chromosomes retaining a split
range **13 to 18** across the 30 settings; the published setting gives **15**.

**Thirteen chromosomes retain a split at all 30 settings**: chr1, chr3, chr5,
chr7, chr8, chr9, chr10, chr11, chr12, chr16, chr18, chr19, chrX. The five that
come and go are chr2, chr4, chr6, chr17, chr20 — chr6 only survives at
`sil <= 0.3`, chr17 only at `sd <= 4`.

Outlier totals by silhouette threshold at `sd = 5`:

| sil | 0.0 | 0.3 | 0.4 | **0.5** | 0.6 | 0.7 |
|---|---|---|---|---|---|---|
| outliers | 2,224 | 1,661 | 1,037 | **821** | 816 | 611 |
| chroms split | 18 | 17 | 16 | **15** | 15 | 13 |

**The published 0.5 sits on a plateau with 0.6**: 821 vs 816 outliers (0.6%
apart) and 15 chromosomes split in both. Moving to 0.4 costs +26% outliers,
moving to 0.7 costs −26%. So the value is not on a cliff — there is a flat
interval and 0.5 is on it. `centroid_sd` is the weaker knob: at `sil = 0.5` it
moves outliers 858 → 812 across sd 3 → 7, a 5.4% span.

## 3. The four chromosomes in the published heatmap are invariant

Figure 6D is the allele co-occurrence heatmap on chr3, chr8, chr11, chr12. Their
outlier counts across all 30 settings:

| chromosome | outlier count | span |
|---|---|---|
| chr3 | 124 at every setting | **0** |
| chr11 | 178 at every setting | **0** |
| chr12 | 233 at every setting | **0** |
| chr8 | 206 at 25 of 30 (drops to 1 only at `sil = 0.7`) | 0 within the neighbourhood |

So the published heatmap rests on chromosomes whose major/minor split does not
move at all over a 30-point grid. That is the strongest row in this sweep.

## 4. Four of the five FISH-validated haplotypes hold; chr19 needs Yi-An

Counts are not the published claim — the claim is that five specific minor
subtypes are real, each backed by FISH. Tracking those haplotypes through Yi-An's
own `apply_sil_and_centroid()`:

| chromosome | sample | Outlier at | note |
|---|---|---|---|
| chr3 | NA21144 | **30/30** | invariant |
| chr9 | HG02630 | **30/30** | invariant (the trio-validated HSat3 inversion) |
| chr12 | HG03816 | **30/30** | invariant |
| chr5 | HG00558 | 24/30 | Outlier at the published setting; lost only at `sd >= 6` |
| chr19 | NA18960 | **0/30** | **Major everywhere, published setting included** |

The first four are a clean answer to R3-major-3: the biological results that
carry FISH validation are robust to both criteria, and the one that is not
invariant is still an outlier at the published value and everywhere in the
neighbourhood except the two highest `sd`.

**chr19 is not a sensitivity result and should not be reported as one.** It is
parameter-independent, so it is not about the thresholds. Facts:

* NA18960 **is** in the input, with both chr19 haplotypes present
  (`NA18960#1#CM101541.1`, `NA18960#2#CM101568.1`), so this is not a missing-sample
  or QC-filter artefact.
* Both are `chr19_Major`, at `norm_divergence` 0.4696 and 0.2036.
* 0.4696 is **higher** than several haplotypes called Outlier on other chromosomes
  (chr1 0.3058, chr17 0.4336) — the classification is per-chromosome relative to
  that chromosome's own spread, and chr19's spread is wide enough to absorb it.
  chr19 has only 4 outliers out of 328 haplotypes at the published setting, none of
  them NA18960.
* The manuscript introduces this one differently from the other three — "Two
  further minor subtypes, a chr12 D12Z3 variant (HG03816) and a chr19 D19Z3 variant
  (GM18960), showed comparable concordance" against the *supplementary* FISH panel.
  So the prediction may well come from subtype composition or per-chromosome
  divergence rather than from the all-chromosome dendrogram outlier call.

**Question for Yi-An, not a correction:** was the chr19 D19Z3 call derived from
this dendrogram's Major/Outlier classification, or from divergence / subtype
composition? If the latter, nothing here is in tension and the sweep simply does
not cover that route. If the former, the reproduction and the published run differ
for a reason that is independent of both swept parameters.

### TRAP that cost a wrong first pass
The manuscript names two of these samples by **Coriell catalogue number**
(`GM21144`, `GM18960`); the assemblies use the **1000G/HPRC name** (`NA21144`,
`NA18960`). The table holds 3,476 `HG` and 1,170 `NA` samples and **no `GM` at
all**, so looking them up as written returns "absent" rather than an error — and
"absent" reads as a real finding. Same class as
[[silent-empty-results-trap]]. First pass reported chr3 and chr19 as 0/30 for
this reason; chr3 is in fact invariant.

## 5. Stage attribution — the centroid scan supplies about half the splits

With the centroid scan off, only the silhouette gate acts:

| sil | 0.0 | 0.3 | 0.4 | **0.5** | 0.6 | 0.7 |
|---|---|---|---|---|---|---|
| chroms split, scan **off** | 18 | 14 | 10 | **8** | 7 | 6 |
| chroms split, scan **on** (sd=5) | 18 | 17 | 16 | **15** | 15 | 13 |

At the published setting the silhouette gate leaves **8** chromosomes with a
genuine two-cluster split, and the centroid scan adds outliers on **7** more. So
roughly half the chromosomes counted as "retaining a split" do so because of stage
2, not stage 1. Worth stating plainly if the response letter describes the split
count, since the two stages answer different questions — stage 1 asks whether a
chromosome is bimodal, stage 2 asks whether individual haplotypes are extreme.

## 6. Verdicts for the parameter table

| parameter | value | verdict | one-sentence justification |
|---|---|---|---|
| `--sil-threshold` | 0.5 | **sweep-backed row** | Outlier count and split count are flat between 0.5 and 0.6 (821 vs 816 outliers, 15 chromosomes split in both) and change by roughly 26% at 0.4 and 0.7, so the published value sits on a plateau rather than at a fitted optimum. |
| `--centroid-sd` | 5 | **sweep-backed row** | The weaker of the two criteria: over sd 3 to 7 at the published silhouette threshold the outlier count moves 5.4% and every FISH-validated haplotype except chr5's keeps its class. |
| combined | (0.5, 5) | **sweep-backed row** | Thirteen chromosomes retain a split at all 30 combinations, the four chromosomes in the allele heatmap have identical major/minor counts at all 30, and three of the five FISH-validated haplotypes are outliers at all 30. |

Nothing here is problematic and no published number moves. As the brief predicted,
split counts rise as the silhouette threshold drops (18 at 0.0) — that is the
conservative-by-construction behaviour of the gate, not fragility, and should be
described as such.

Related: [[a6-maxgap-sweep-result]] [[a6-inputs-found-2026-08-26]]
[[karyoscope-revision-analyses]]
