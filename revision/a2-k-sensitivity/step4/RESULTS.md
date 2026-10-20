# A2 Step 4 — the k-robustness questions Step 3 left open

R3-major-3 asks four things. Step 3 (`../README.md`) answered the Robertsonian/SST1
breakpoint: invariant across *k* = 21–61. This directory answers the other three —
the **novel fraction**, **performance in divergent regions**, and the two remaining
biological applications, **D4Z4 array reconstruction** and **centromere subtype
clustering**.

**They do not behave the same way, and two of them are negative results.**

| analysis | robust to *k*? | one line |
|---|---|---|
| novel fraction | n/a (a measurement) | rises ~0.11 pp per unit of *k*; median novel run == *k* |
| divergence simulation | yes, smoothly | *k* = 31 never more than 1.7 points off the best at any tested divergence |
| **D4Z4 census** | **no** | *k* = 51 returns **zero arrays**; *k* = 21 halves the count. **Salvageable** by scaling the caller's two base-pair thresholds with *k* |
| **centromere clustering** | **no** | only 67.9% of haplotypes keep their call between *k* = 21 and *k* = 51; 2 of 4 FISH-validated subtypes lost off *k* = 31 |

Full write-up with the framing recommended for the response letter:
`../../docs/write-ups/results-handoff-6-2026-09-03.md`.

---

## The database trap — read this before reproducing anything

**Both biological analyses consume *v1* (`KS_human_CHM13`) region labels, and the
v1-vs-v2 difference is larger than anything *k* does.** The six-*k* family built for
Step 3 (`HKS_human_CHM13_v2_k*`) is v2 and cannot be used for either:

* **D4Z4**: v1's region source labels the 4q35/10q26 arrays plain `q_arm`, so each
  3.3 kb unit picks up `bsat` only where its k-mers are shared with beta-satellites
  elsewhere — and that alternation of `bsat_specific` / `q_arm_specific` 100 bp bins
  is exactly what Floris's caller counts. **v2 labels the whole 117 kb array as one
  `bSat` interval**, so every array collapses to one unit: 28,580 → 1,378 at *k* = 31.
* **centromere**: at *k* = 31 the v2 family gives 1,195 minor haplotypes against the
  published 821 — the vocabularies differ, so the feature matrix differs.

Fix: the v1 region set was regenerated exactly and rebuilt with HKS at each *k*
(`scripts/v1-region-rebuild/`). Recipe:

1. Input = `KaryoScope-archive:resources/KS_human_CHM13/old.CHM13.censat.processed.bed`
   **minus rows labelled `rDNA`** (archive commit `cbc091d` removed them; the dropped
   9,928,083 bp become interior gaps the generator labels `noncentromeric`).
2. Generator = `KaryoScope-archive:workflow/scripts/get_region_featureset.py` **@ `0d1cb76`**
   (2025-08-25 — the version that built v1). *Not vendored here*, per repo policy;
   it is Snakemake-driven, so `run_v1_generator.py` supplies the stub object.
3. Verified: leaf set identical to the installed v1 `hierarchy.tsv`; priority file
   identical to v1's `KS_human_CHM13.priority.txt` region rows; tiles the genome with
   no gaps (2271 rows, 3,117,292,070 bp = the whole genome incl. chrM as `exclude`).
   The priority file is committed as `results/CHM13_reg.v1.priority.txt`
   (sha256 `bd78bbeb262088cb…`); **the BED is not — this repo's `.gitignore` excludes `*.bed`** —
   but it is deterministic, sha256 `5acf54ac91b04f22…`. Regenerate exactly:

   ```bash
   ARCH=<KaryoScope-archive>
   awk -F'\t' '$4!="rDNA"' $ARCH/resources/KS_human_CHM13/old.CHM13.censat.processed.bed \
       > CHM13.censat.processed.v1.bed
   git -C $ARCH show 0d1cb76:workflow/scripts/get_region_featureset.py > gen_v1.py
   python3 scripts/v1-region-rebuild/run_v1_generator.py \
       $ARCH/resources/KS_human_CHM13/CHM13/KS_human_CHM13.CHM13.fasta.gz.fai \
       CHM13.censat.processed.v1.bed CHM13_reg.v1.bed CHM13_reg.v1.priority.txt gen_v1.py
   ```
4. Build with `build.k{21,31,51}.yaml`. **Use `build: external_memory:`** — 20 GB peak
   instead of ~75 GB, and the jobs schedule immediately at 40 GB. Measured on Ice Lake,
   16 threads: 15 min (*k* = 21), 17.5 min (*k* = 31), 36.5 min (*k* = 51); peak RSS
   19.7 / 19.8 / 20.5 GB.

`bin_features_feb2026.py` in `scripts/v1-region-rebuild/` is a **pinned copy** of
`KaryoScope-archive:workflow/scripts/bin_features.py` @ `18dccc6` — the version that
produced the published inputs. The archive has since gained a ">50% novel" bin rule
(`ea5c2fe`) and the v2 vocabulary (`f15c29b`), neither of which was in effect then.

**Both controls came out exact**, which is what licenses the *k* comparisons:
the installed v1 KMC database on the same windows reproduces the published D4Z4
census on **582/582 haplotypes in every column**; the v1-like HKS *k* = 31 rebuild
reproduces that **bit-for-bit**; and the v1-like *k* = 31 centromere clustering
reproduces the published run on **4,551/4,551 haplotypes**, same call *and* same
cluster label.

## Job 1 — novel fraction vs *k* (`results/novel_by_k.tsv`)

No re-annotation; reads the Step 2 BEDs. Method mirrors A4 novel-v3 (N runs masked
with the *k*−1 left extension, chrM excluded). Anchor reproduced exactly:
HG002 hap1 at *k* = 31 gives 75,396,753 bp = 2.5643% of 2,940,260,101 evaluated
positions, 86.90% in runs <= *k*, 0.60% in runs > 1 kb. Pre-smoothing novel is
identical across all six feature sets at every *k*.

HG002 hap1: 1.45 / 2.56 / 3.70 / 4.81 / 5.87% at *k* = 21/31/41/51/61.
CHM13 self-annotation: **0 novel bp at every *k***. The median novel run is exactly
*k* for *k* >= 31 — the signature of single-site divergence, not absent sequence.

## Job 2 — divergence at *k* = 21/31/51 (`results/divergence_by_k.tsv`)

A5 design reused exactly; the four mutated genomes were regenerated and verified
byte-identical to A5's. **Control: the *k* = 31 rows are identical to A5's shipped-database
rows in all 8 cells.** Measured novel tracks 1−(1−r)^k at every rate, always slightly
below it. Post-smoothing agreement is much flatter than pre-smoothing, and the best
*k* depends on divergence: *k* = 51 at 0.1%, *k* = 31 at 1% (95.5%), and *k* = 51
collapses at 5% (66.4% vs 83.8%).

## Job 3 — D4Z4 (`results/d4z4_by_k.tsv`, `d4z4_by_k_summary.md`)

| *k* | total units | main-array units | arrays | HMM Pearson *r* (n = 231) |
|---|---|---|---|---|
| 21 | 13,620 | 10,193 | 624 | 0.9424 |
| **31** | **28,580** | **14,433** | **629** | **0.9977** |
| 51 | 33,332 | **0** | **0** | undefined |
| published | 28,580 | 14,433 | 629 | 0.9977 |

**Mechanism.** *k* = 21 halves the count by fusing fragments (78% of inter-fragment
gaps fall below the caller's 400 bp merge threshold, vs 66% at *k* = 31). *k* = 51
fails on the **size gate, not the exon gate**: median merged BSR cluster falls
2,700 → 800 → 200 bp, so at *k* = 51 only 0.8% of clusters clear the hard-coded
`MIN_ARRAY_ELEMENT_SIZE` = 1,500 and nothing is ever promoted to an array element.
The exon gate passes 84–89% of clusters at every *k*.

**Salvage: yes** (`results/d4z4_salvage_sweep.tsv`). Scored against the HMM counts,
which are external and identical at every *k* — not against the published totals,
which would be circular.

| | published (*k*=31) | best *k*=21 | best *k*=51 |
|---|---|---|---|
| setting | merge 400, min 1500 | merge 300, min 100 | merge 800, min 100 |
| main-array units | 14,433 | 14,244 | 14,376 |
| arrays | 629 | 630 | 630 |
| HMM *r* | 0.9977 | 0.9952 | 0.9970 |
| total units | 28,580 | 18,156 | 29,848 |

The best merge threshold **scales with *k*** (300 / 400 / 800), as the mechanism
predicts for a fixed base-pair gap applied to a signal whose granularity *k* sets.
What does not recover at *k* = 21 is `repeat_count_total` (the published 28,580
headline), which counts every BSR cluster genome-wide.

**Fig. S8B is now fully reproducible.** The HMM count table
(`HPRC-D4Z4:Hailey-notes/HPRC_chr4_counts.xlsx`, 534 contigs) reproduces the published
*r* = 0.998 at *n* = 231 exactly. An independent re-run of the Methods HMM recipe on
the same windows (`scripts/job3-d4z4/run_hmm_counts.sh`) agrees with that table at
*r* = 0.9991.

## Job 4 — centromere clustering (`results/centromere_by_k.tsv`)

| *k* | minor | chromosomes with a split | concordance vs *k*=31 |
|---|---|---|---|
| 21 | 395 | 16 | 76.4% |
| **31** | **821** | **15** | — (= published, 4,551/4,551) |
| 51 | 1,266 | 17 | 85.9% |

**Only 67.9% of haplotypes keep their call between *k* = 21 and *k* = 51.** The four
Fig. 6D chromosomes — invariant across A6's 30-point gate grid — collapse at *k* = 21
(chr3 124 → 2 minor, chr8 206 → 1). **Two of the four FISH-validated subtypes are lost
off *k* = 31**: chr3 NA21144 (both haplotypes major at *k* = 21) and chr9 HG02630
(major at *k* = 21 *and* *k* = 51). chr5 HG00558 and chr12 HG03816 hold at every *k*.

**Mechanism — loss of leaf resolution, not fragmentation.** The bp composition is
nearly *k*-invariant; the *label vocabulary* is not. At low *k* the alpha-satellite
subtypes share k-mers, so the hierarchy-aware assignment can only report their common
ancestor `asat`:

| *k* | leaf-labelled bp | internal-node bp | `asat` alone | `active` |
|---|---|---|---|---|
| 21 | 81.3% | 18.7% | 14.1% | 3.4% |
| 31 | 92.2% | 7.7% | 0.2% | 11.6% |
| 51 | 97.0% | 2.6% | 0.01% | 16.3% |

NA21144's chr3 haplotype over 90.3–91.5 Mb:

```
k=21:  monomeric  asat  divergent  asat  inactive  asat  ct
k=31:  monomeric        divergent        inactive        ct
```

Two effects compound: the subtype distinctions the clustering separates on disappear,
**and** the unresolved stretches inject junk transitions. **KaryoScope is behaving
correctly** — 21-mers genuinely do not separate active from inactive HOR — so at
*k* = 21 this is information loss, not a representation artefact.

**Salvage: partial, and not worth adopting** (`results/centromere_representation_salvage.tsv`).
Five feature representations x 3 *k*, gated exactly as published:

| representation | 21 vs 31 | 51 vs 31 | 21 vs 51 | *k*=31 vs published | FISH minor (of 12) |
|---|---|---|---|---|---|
| **published** (edges, blockweight) | 76.4% | 85.9% | **67.9%** | **100%** | 9 |
| length-weighted edges | 82.9% | 81.9% | 74.0% | 86.3% | 7 |
| abundance only | 81.1% | 78.4% | 79.7% | 81.4% | **10** |
| abundance only, length-weighted | **83.6%** | **87.5%** | 81.2% | 79.1% | 7 |
| edges + abundance | 80.0% | 80.8% | **85.1%** | 81.9% | **10** |

Cross-*k* concordance at the extremes improves from 67.9% to ~85%, but every
alternative changes the published *k* = 31 calls by 14–21% of haplotypes, and **the
representations that most improve *k*-stability are the ones that lose FISH-validated
subtypes** — they buy agreement by flattening real structure. One robust observation:
the number of chromosomes retaining a split is **16 / 15 / 17 in all five
representations**, so the chromosome-level structure is representation-independent.

---

## Two method errors of mine, caught and fixed

1. **Centromere gating.** The first pass gave 833 minor against the published 821.
   `KS_allchr_barplot.py` defaults to `--exclude-features novel` and drops
   `EXCLUDED_CHROMS` **before** gating; I passed no exclusion and dropped them after.
   Corrected -> exact reproduction. Always mirror that script's own defaults.
2. **The D4Z4 sweep's merge dimension never varied.** See
   `../../a6-parameter-sweeps/d4z4/RESULTS.md` — this also forces a correction to A6.

## Inputs not stored here

* `release2_contigs.tsv` (4.3 MB), the contig -> length + FASTA index. Regenerate:
  ```
  for f in $KS_REFERENCES/pangenome/release2/*.fa.gz.fai; do
    awk -v F="${f%.fai}" 'BEGIN{OFS="\t"}{print $1,$2,F}' "$f"; done > release2_contigs.tsv
  ```
* The window FASTAs (291 Mb for D4Z4, 72.3 Gb per *k* for centromeres) — extracted
  with `samtools faidx` from the coordinates in `results/{d4z4,cen}_windows.tsv`.
* `HPRC_chr4_counts.xlsx` — Hailey Loucks's HMM counts, from `barthel-lab/HPRC-D4Z4`.
* The v1-like databases themselves, under
  `$KS_WORK/ks_build_sources/v1like_ksweep/db`.

See `../../docs/INPUTS.md`.
