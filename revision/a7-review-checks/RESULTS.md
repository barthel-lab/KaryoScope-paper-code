# A7 — five checks raised by the 2026-09-07 independent review

Brief: [`docs/briefs/a7-review-checks-brief.md`](../docs/briefs/a7-review-checks-brief.md). Scripts as run in
`scripts/` (SLURM directives included; `job2b_noprio.build.yaml` is the gene-only, no-priority build spec),
result tables in `results/`, per-job provenance (versions, nodes, wall-clock, peak RSS, SLURM ids) in
`provenance/`. Cluster working copy: `$KS_WORK/a7_review_checks/`.


Cluster session, 2026-09-07.
All numbers below are measured; nothing is quoted from the manuscript except where marked "published".
Software: karyoscope 2.3.0, hks 0.4.0, HKS_human_CHM13_v2 (shipped registry), HKS_arabidopsis_ColCEN 1.2.0.

Deliverables (all in `results/`):

| file | job |
|---|---|
| `ey15_2_novel_by_stage.tsv` (+ `ey15_2_labels_by_stage.tsv`, `ey15_2_novel_fate.tsv`, `ey15_2_novel_sanity.tsv`, `job1_provenance.txt`) | 1a |
| `ey15_2_chromosome_by_sequence.tsv` | 1b |
| `hg002_gene_labels_by_stage.tsv` (+ `chm13_gene_labels_by_stage.tsv` control, `hg002_gene_labels_by_stage_NOPRIORITY.tsv`) | 2 |
| `chm13_gene_db_composition.tsv` | 2 |
| `rc_symmetry_hks.tsv`, `rc_symmetry_kmc.tsv` (+ `_label_pairs.tsv`, `_run_hist.tsv`) | 3 |
| `acro_multilabel_flanks.tsv` (+ `context_chr13_9202310_{pre,post}.tsv`) | 4 |
| `d4z4_multiarray_contigs.tsv` | 5 |

## Job 1 — Ey15-2: the pre-smoothing novel fraction is 18.66 %, not 105 bp

No re-annotation: A3 job 35378568 kept both stages (`a3_crossgenome/annot_ey15_2/*.presmoothed.bed.gz`).
Denominator = tiled k-mer positions, N-excluded per A4 (one 100-bp N run on OX298799.1 at 17,215,895, extended 30 bp
left; 130 bp masked). Col-CEN nuclear-only self-annotation is the control (no N).

**1a — `ey15_2_novel_by_stage.tsv`** (denominator 136,135,091 bp Ey15-2; 131,559,526 bp Col-CEN):

| set | stage | label-specific % | multilabel % | novel % (bp) |
|---|---|---|---|---|
| chromosome | presmoothed | 68.59 | 12.75 | **18.66 (25,402,658)** |
| chromosome | smoothed | 89.65 | 10.35 | 0.00008 (105) |
| region | presmoothed | 72.71 | 8.63 | 18.66 |
| region | smoothed | 96.28 | 3.16 | 0.558 (760,084) |
| repeat | presmoothed | 69.98 | 11.37 | 18.66 |
| repeat | smoothed | 95.46 | 2.70 | 1.834 (2,497,066) |
| gene | presmoothed | 69.88 | 11.46 | 18.66 |
| gene | smoothed | 93.48 | 2.98 | 3.544 (4,824,275) |
| Col-CEN (all four sets) | presmoothed | 88.1–91.6 | 8.4–11.9 | **0** |
| Col-CEN (all four sets) | smoothed | 99.3–99.996 | 0.004–0.66 | 0 |

* **Sanity passed to the base**: the smoothed novel column on the unmasked BEDs reproduces `label_composition.tsv`
  exactly — 105 / 760,084 / 2,497,066 / 4,824,275 (`ey15_2_novel_sanity.tsv`). The brief's "2,496,xxx" for repeat
  is a typo; the published table says 2,497,066. The 105 bp is not inside the N-mask, so masking changes nothing.
* **The representation figure is 25.4 Mb = 18.66 % of Ey15-2's k-mer positions**, identical across the four sets
  (as it must be: a k-mer is either in the base index or not). 105 bp is the imputation residue of the chromosome
  set. PR #8's "essentially none of the query is unrepresented" is wrong by that definition; the correct
  statement is that 18.66 % of Ey15-2 k-mers are absent from the Col-CEN index (vs 0 % for Col-CEN itself), and
  smoothing imputes a label for all but 105 bp of them.
* **Where the 25.4 Mb goes** (`ey15_2_novel_fate.tsv`): chromosome set — 77.2 % to a specific chromosome, 22.8 %
  to `nuclear`, 0 % left novel; region 87.0 % leaf / 10.0 % internal / 3.0 % novel; repeat 83.7 / 6.5 / 9.8;
  gene 81.0 / 0 / 19.0 (the gene set has no internal node below the root, so unresolvable bases stay novel).
* Multilabel (`nuclear`) is 12.75 % pre and still 10.35 % post for Ey15-2, against 11.9 % → 0.004 % in Col-CEN:
  in the self-annotation the cross-chromosome k-mers are resolved by their flanks; in Ey15-2 the flanks disagree.

**1b — `ey15_2_chromosome_by_sequence.tsv`.** Chromosome identity was read from the assembly record's own FASTA
headers (`OX298797.1 ... chromosome: 1` … `OX298801.1 ... chromosome: 5`) — the accessions do run 1→5 in order, but
that was checked, not assumed. Col-CEN and Ey15-2 both use TAIR Chr1–Chr5; no numbering difference.

| sequence (ENA chr) | stage | own chromosome % | other chromosome % | nuclear % | novel % |
|---|---|---|---|---|---|
| OX298797.1 (Chr1) | pre / post | 68.6 / **87.7** | 0.86 / **2.77** | 11.8 / 9.5 | 18.7 / 0 |
| OX298798.1 (Chr2) | pre / post | 65.9 / **82.5** | 1.16 / **5.35** | 13.7 / 12.1 | 19.3 / 0 |
| OX298799.1 (Chr3) | pre / post | 66.6 / **85.5** | 1.16 / **3.82** | 13.7 / 10.6 | 18.6 / 0 |
| OX298800.1 (Chr4) | pre / post | 66.4 / **86.2** | 1.16 / **3.77** | 13.1 / 10.1 | 19.3 / 0 |
| OX298801.1 (Chr5) | pre / post | 68.9 / **86.8** | 1.25 / **3.27** | 12.1 / 9.9 | 17.8 / 0 |
| Col-CEN Chr1–5 (control) | pre / post | 85.4–90.2 / 99.99–100 | 0 / 0 | 9.8–14.6 / ≤0.01 | 0 / 0 |

* "Correct but less specific" is **partly** right: 9.5–12.1 % of each sequence ends as `nuclear` (correct, less
  specific). But **2.8–5.4 % of each Ey15-2 chromosome is labelled as a different chromosome after smoothing**, up
  from 0.9–1.3 % before it — smoothing extends the cross-chromosome k-mer calls, it does not suppress them. In the
  Col-CEN control the wrong-chromosome fraction is 0 at both stages. If PR #8 keeps the hierarchy sentence, it
  needs this number next to it. (Whether some of that 3–5 % reflects real Ey15-2 vs Col-0 structural differences
  was not tested here.)

## Job 2 — HG002 hap1 gene set: the numbers are right, the labels are permuted

`hg002_gene_labels_by_stage.tsv`, from the A4 BEDs (HKS_human_CHM13_v2, job 35378688), A4 N-mask, denominator
**2,940,260,101**. Control: CHM13 self-annotation (`chm13_gene_labels_by_stage.tsv`, no mask, chrM excluded,
denominator 3,117,274,781).

| label | HG002 hap1 pre % | HG002 hap1 smoothed % (bp) | CHM13 self, both stages % | build input BED % (reference bp) | index k=31 node count % |
|---|---|---|---|---|---|
| exon | 11.38 | **11.47** (337,237,167) | 11.56 | **4.77** | 5.55 |
| intron | 48.43 | **49.54** (1,456,574,163) | 49.05 | 48.77 | 52.68 |
| intergenic | 37.62 | **38.76** (1,139,527,175) | 39.39 | 46.46 | 41.77 |
| categorized | 0 | 0 | 0 | 0 | 0 |
| novel | 2.56 | **0.2354** (6,921,596) | 0 | — | — |

* **Results L195 is a label permutation.** "exon / intron / intergenic = 49.4 / 38.8 / 11.6 %" — the measured
  smoothed values are **intron 49.5 / intergenic 38.8 / exon 11.5 %**. The three numbers are right to ≲0.15
  points; they were attached to the wrong labels (handoff #2 §7 carried the same order, so the error is mine, not
  the manuscript session's). Exon is 11.5 %, not 49.4 %. The 0.2354 % novel post reproduces exactly.
* **The propagation effect is real but 2.4×, not 10×.** Exons are 4.77 % of T2T-CHM13v2.0 bases in the build
  input (`chm13_gene.bed`, sha256 03ae30a1…, 612,909 rows, tiles the genome incl. chrM), yet 11.56 % of CHM13's own
  k-mer positions are labelled `exon` after priority resolution — every genomic copy of a k-mer that occurs in any
  exon is labelled exon. Intergenic loses correspondingly (46.5 → 39.4 %); intron is nearly neutral (48.8 → 49.1 %).
  The Methods sentence should be written around 4.8 % → 11.5 %, not around 49 %.
* Index composition (`chm13_gene_db_composition.tsv`, from A4's `hks node-stats` at k=31): exon 278.6 M / intron
  2,647.2 M / intergenic 2,098.9 M nodes, `categorized` 0 — priority resolution leaves no tied k-mer. The index
  stores both strands (hks default; no `--forward-only`), so node counts are ~2× distinct canonical k-mers; the
  fractions are what matter. The CHM13 self-annotation is bit-identical pre and post smoothing (no novel, no
  multilabel — nothing to smooth).
* **Without the priority file** (`hg002_gene_labels_by_stage_NOPRIORITY.tsv`; gene-only build from the same input
  BED, 15 min / 19.7 GB with the spec-level external-memory block, then HG002 hap1 annotate): presmoothed
  `categorized` = **396,600,461 bp (13.49 %)**, exon 126.8 Mb (4.31 %), intron 1,235.2 Mb (42.01 %), intergenic
  1,106,258,434 bp (37.62 %) — intergenic is **bit-identical** to the priority build's presmoothed value, as it
  must be (lowest priority: it can only ever be the sole claimant). With priorities the 396.6 Mb splits ~208 Mb to
  exon and ~189 Mb to intron. Smoothing then resolves most of it (categorized 0.77 % post; exon 4.69 %).
  **The recipe's "786 Mb on categorized, exon 258 Mb vs 679 Mb" does not reproduce on hap1**: the measured hap1
  values are 396.6 / 126.8 / 334.5 Mb, each 1.98–2.03× smaller — the docs numbers are consistent with the combined
  diploid `hg002v1.1.fasta.gz` (both haplotypes), which was not rerun here. Quote either the diploid with that
  qualifier, or the hap1 figures above.

## Job 3 — reverse-complement symmetry: exact on the HKS path

`rc_symmetry_hks.tsv`. Queries: HG002 hap1 `chr21_hap1` (44,312,598 bp, acrocentric), `chr20_hap1` (67,035,379 bp,
ordinary autosome), and the rob(14;21) fusion scaffold of one Robertsonian sample. [Snapshot note: the fusion
scaffold's name and length, its rows in the `rc_symmetry_*` tables and the Job 3 scripts are withheld pending the
dbGaP ruling on the ROB assemblies; the HG002 rows are complete.] Reverse-complemented with `seqtk seq -r`; verified that
rc(rc(x)) = x AND that the RC head equals the RC of the forward tail (i.e. the file really was reversed). Both
orientations annotated in one `karyoscope annotate` run, HKS_human_CHM13_v2, k = 31, all six feature sets,
`hks smooth --max-gap 1000` (annotate's default; **on the HKS backend the smoother is `hks smooth`, src/smooth.rs,
not `smooth.py`** — the KMC run below exercises `smooth.py`). RC intervals mapped back as [L'−e, L'−s) with
L' = L−k+1; the presmoothed stage is the control for the mapping.

| backend | sequences × sets × stages | bp compared per stage | bp differing | differing runs |
|---|---|---|---|---|
| HKS (hks smooth) | 3 × 6 × 2 = 36 | 244,285,671 | **0 in every comparison** | 0 |
| KMC (smooth.py) | 3 × 6 × 2 = 36 | 244,285,671 | **0 in every comparison** | 0 |

* Presmoothed identical (control passed — both strands are in the index). **Smoothed identical to the base in all
  18 comparisons**: the labelling of a sequence and of its reverse complement are not just "equivalent", they are
  the same after coordinate reflection, on the acrocentric chromosome, the ordinary autosome and the fusion
  scaffold. The RC output BEDs differ byte-wise from the forward ones (different md5, same row count), so this is
  not a trivial identity.
* **The KMC database (`KS_human_CHM13_v2`), whose smoothing runs in `karyoscope/core/smooth.py` — the "Advance i"
  scan the brief reasons about — gives the same result: 0 differing bases in all 36 comparisons.** So the
  asymmetric advance rule does not produce an asymmetric labelling on these three sequences; the Methods sentence
  can be stated as "produces identical labellings" (measured), not merely "designed so that". Bonus: the KMC and
  HKS presmoothed BEDs are byte-identical (same md5 on chr21 region), re-confirming A2's KMC ≡ HKS.
* Runtime: six-input HKS annotate 3 min 19 s / 10.4 GB peak; KMC 24 min 04 s / 18.9 GB.

## Job 4 — long acrocentric multilabel runs: the flank rule holds, with one refinement

`acro_multilabel_flanks.tsv`, from the A4 CHM13 and HG002 hap1 acrocentric BEDs (HKS_human_CHM13_v2). Hierarchy:
categorized → {acrocentric → {array → {DJ, PJ, rDNA}, SST1, PHR}, nonacrocentric}. A run = a maximal stretch of
consecutive presmoothed intervals that are not leaves (internal labels or `novel`); reported when its internal-label
bp > 1 kb. Every such run in both assemblies is labelled **`array`** (no `acrocentric` or `categorized` run exceeds
1 kb) and every one lies on chr13/14/15/21/22, inside the rDNA arrays.

| assembly | runs > 1 kb | multilabel bp | longest | immediate-flank rule consistent | effective-flank rule consistent |
|---|---|---|---|---|---|
| CHM13 | 427 | 1,184,009 | 5,656 | 417 / 427 | **427 / 427** |
| HG002 hap1 | 68 | 170,115 | 5,493 | 66 / 68 | 61 / 68 |

* **Immediate flanks** (the leaf intervals directly adjacent to the run): CHM13 414 runs flanked rDNA/rDNA → all
  reassigned to rDNA; 7 PJ/PJ → PJ; 3 rDNA/PJ (LCA = `array` = run label) → unchanged. HG002: 37 rDNA/rDNA and
  3 PJ/PJ reassigned; 26 with LCA not below `array` (rDNA/PJ, PHR/*, nonacrocentric/*) unchanged. **12 runs
  contradict the immediate-flank rule** — all reassigned to rDNA although one or both adjacent leaves were a
  `PHR` interval of 10–47 bp (single k-mer blips). Example CHM13 chr13:9,202,310–9,204,039 (`context_chr13_*`):
  `rDNA … array(470) PHR(31) array(1729) PHR(24) array(769) rDNA(31) …` → everything but the PHR blips becomes rDNA.
* **Refinement, mirroring the window scan**: an unrelated-branch interval (PHR, SST1, nonacrocentric) inside
  max_gap is skipped, not treated as a flank; the effective flanks are the nearest leaves in the run label's own
  subtree, or an interval at least as general as the run label, which caps the window. Under that rule CHM13 is
  427/427 consistent (417 rDNA/rDNA and 7 PJ/PJ reassigned; 3 PJ/rDNA unchanged). Skipped unrelated spans were
  ≤ 78 bp — max_gap never came into play.
* **The 7 HG002 exceptions** (all "expected rDNA, stayed `array`", chr15 ×6, chr22 ×1) share one signature: a
  `novel` run (3–27 bp) sits at the run's edge next to an off-branch blip (`PHR` 2–28 bp or `nonacrocentric`
  11 bp). That novel run is imputed to LCA(PHR, array) = `acrocentric` (or, with nonacrocentric, stays at the
  root), and that less-specific interval then bounds every window that reaches it, so the neighbouring `array`
  runs never see their rDNA flanks. E.g. chr15_hap1 5,850,599: `PHR(28) novel(3) array(15) rDNA(6)` → post
  `PHR(28) acrocentric(3) array(15) rDNA…`, and the 3,243-bp and 2,392-bp `array` runs to its left stay `array`.
  It is an order-of-operations effect of the fixed-point pass, not a max_gap effect. CHM13 has no such case.
* So the Methods sentence should read: long multilabel runs are reassigned when their nearest same-branch flanks
  share an ancestor more specific than the run's label (rDNA/rDNA → rDNA), and survive when they do not (a PJ
  flank on one side and rDNA on the other: LCA = `array`, the run's own label — i.e. the run genuinely sits
  between two distinct arrays) or when an adjacent novel blip was first imputed to a less specific label.
  `max_gap` bounds skipped unrelated-branch spans, not run length, and never bit here.

## Job 5 — the 25 multi-array D4Z4 haplotypes are same-contig by construction

`d4z4_multiarray_contigs.tsv`: 80 array rows over 28 haplotypes — the **25 multi-array** haplotypes (20 canonical-only,
5 mixed canonical+degraded) plus the **3 degraded-only single-array** haplotypes, so the categories' overlap is
visible (`class` column). Source: the frozen k = 31 caller table
`a2_step4/job3_d4z4/base_published/results/d4z4_array_analysis.tsv`, `cmp`-identical to the v1-like k = 31 rerun.

* **`all_arrays_same_contig` is TRUE for all 25, trivially.** The caller receives one 500-kb q-terminal window of
  one contig per haplotype (`common/d4z4_windows.tsv`, 582 rows, one contig each; the caller filters its
  annotations on `c == hap_name`, where the name embeds the contig accession), so no haplotype can have arrays on
  two contigs. **The response-letter sentence is a property of the input design, not an observation** — it should
  say that only T2T q-arm contigs entered the census, so a contig break cannot split an array. The stronger
  observable check is `is_truncated`, which is False for all 28.
* Per-array columns: window-frame and contig-frame coordinates (`array_start_contig` = window start + x), element
  count, type, terminal label (pLAM 5.2–5.8 kb / ct_terminal / none), pLAM end, and `gap_from_previous_array_bp`.
  The 52 gaps between consecutive arrays fall in three groups: 7,500–8,300 bp (27 gaps; 25 of them exactly 7,500 or
  7,600), 15,800–20,200 bp (14) and 37,900–94,500 bp (11). The three chr4 mixed haplotypes have a degraded array
  42.6–72.6 kb proximal to a ct-terminal canonical array.

## Surprises, in order of consequence

1. **Job 2**: the manuscript's exon/intron/intergenic percentages are permuted — exon is 11.5 %, not 49.4 %.
   The propagation effect to describe is 4.8 % → 11.5 %.
2. **Job 1**: Ey15-2 has 18.66 % unrepresented k-mers; smoothing imputes 77 % of them to a chromosome, and 2.8–5.4 %
   of each chromosome ends up labelled as another chromosome (0 in the self-control).
3. **Job 3**: smoothing is exactly RC-symmetric on BOTH backends (0 bp in 244 Mb × 6 sets × 2 stages, hks smooth and smooth.py alike).
4. **Job 4**: the flank rule is right once "flank" means nearest same-branch label; the 7 exceptions are novel
   blips imputed to a general label first.
5. **Job 5**: same-contig is by construction; say so.

## Traps hit

* `karyoscope build --external-memory DIR` (CLI flag) did **not** run the k-mer sort in external memory — the
  build hit 50 GB in 55 s and was OOM-killed at 48 GB (`hks build-base` ran with `--mem-gigas 8` but in-memory
  sorting). The spec-level `build: {external_memory: DIR, mem_gigas: 16}` block (as used for the v1-like rebuilds)
  does: same build then peaked at 19.7 GB in 15 min. **Verified cause**: `commands/build.py::_build_spec` returns
  `BuildSpec.from_yaml(spec_path)` as soon as `--spec` is given, so `--external-memory`, `--mem-gigas`, `--threads`
  and `--forward-only` are silently ignored in spec mode (the recorded `hks build-base` command had no
  `--external-memory` and `-t 4 --mem-gigas 8`). The published recipe (`docs/recipes/human-chm13v2.md`, and my
  `recipe_repro/repro.sbatch`) uses exactly this `--spec ... --external-memory` combination and only succeeded
  because it asked for 160 GB. Bug to file on KaryoScope: either merge the flags into the YAML spec or refuse the
  combination.
* `hks node-stats` counts both strands; do not read its counts as distinct canonical k-mers.
