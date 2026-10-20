# Ey15-2 panel label composition — 2026-08-24

Table: `label_composition.tsv` (run, feature_set, label, bp, pct).
Script: `../label_composition.sh`.

Computed from the **same scaffolded+smoothed BEDs the karyotype panels were
rendered from**, so the fractions describe what is actually drawn. Annotations
tile the input, so per-label bp sums exactly to the genome length — verified:

| run | tiled total | note |
|---|---|---|
| Ey15-2 | 136,135,221 bp | nuclear only |
| Col-CEN | 131,559,526 bp | **nuclear only** — organelles are absent from the scaffolded BEDs, so this is a clean nuclear-vs-nuclear comparison despite Col-CEN having been annotated from the full assembly |

Percentages are of each genome's own length; the genomes differ by 4.6 Mb.

## CEN180 — the number asked for

**CEN180 lives in the `region` feature set** (not `repeat`), which is the set the
`centromere.region` and `subtelomere.region` panels render.

| label | Ey15-2 bp | % | Col-CEN bp | % |
|---|---|---|---|---|
| **CEN180** | **13,689,717** | **10.056** | **11,744,336** | **8.927** |
| q_arm | 68,090,453 | 50.017 | 69,763,803 | 53.028 |
| p_arm | 47,249,904 | 34.708 | 47,201,038 | 35.878 |
| arm (generic) | 3,640,446 | 2.674 | 3,273 | 0.003 |
| cen_gap | 2,039,905 | 1.498 | 2,820,194 | 2.144 |
| novel | 760,084 | 0.558 | 0 | 0.000 |
| categorized | 530,376 | 0.390 | 1,679 | 0.001 |
| centromeric | 134,336 | 0.099 | 25,203 | 0.019 |

Ey15-2 carries **~1.9 Mb more CEN180** than Col-CEN (10.06% vs 8.93%),
consistent with the accession being chosen for its resolved, divergent
centromeres.

## The pattern to caption, in one sentence

**Cross-genome, mass moves from specific leaf labels to their generic PARENT
labels — not to `novel`.** Compare the generic label in each set:

| set | generic label | Ey15-2 | Col-CEN | ratio |
|---|---|---|---|---|
| region | `arm` | 2.674% | 0.003% | ~890x |
| chromosome | `nuclear` | 10.349% | 0.004% | ~2600x |
| gene | `categorized` | 2.976% | 0.253% | ~12x |
| repeat | `Transposable_Element` | 1.283% | 0.212% | ~6x |

**This is the smoothing hierarchy working as designed**, and it is the honest
cost of annotating an unseen accession: the sequence is represented, it just
cannot be resolved to the finest label.

### The chromosome set is the sharp example — do not caption it as "0% novel"

| label | Ey15-2 % | Col-CEN % |
|---|---|---|
| Chr1 | 22.742 | 24.734 |
| Chr5 | 20.314 | 22.409 |
| Chr3 | 17.546 | 19.568 |
| Chr4 | 14.863 | 16.400 |
| Chr2 | 14.186 | 16.886 |
| **nuclear** (unresolved) | **10.349** | 0.004 |
| novel | 0.000 (105 bp) | 0.000 |

`novel` really is ~zero (105 bp of 136 Mb), but **10.3% of Ey15-2 resolves only
to `nuclear`, not to a chromosome.** Saying "chromosome assignment transfers
with 0% novel" is true and misleading; the accurate statement is *"essentially
nothing is unrepresented (0.00008% novel), and 89.7% resolves to a specific
chromosome, with the remaining 10.3% falling to the parent `nuclear` node."*

## Other sets, for completeness

**repeat** (Ey15-2 % / Col-CEN %): nonrepeat 64.521 / 67.884, DNA_other 11.632 /
9.759, Helitron 6.914 / 7.930, Gypsy 4.944 / 5.308, Mutator 3.034 / 3.491,
**novel 1.834 / 0**, Transposable_Element 1.283 / 0.212, Copia 1.131 / 1.413,
CACTA 0.960 / 1.067, L1 0.786 / 0.985, Class_II_DNA 0.684 / 0.223,
LTR_other 0.541 / 0.238, categorized 0.436 / 0.182, then 13 labels below 0.35%.

**gene**: intergenic 41.223 / 42.918, exon 34.455 / 36.318, intron 17.802 /
20.511, **novel 3.544 / 0**, categorized 2.976 / 0.253.

### Bonus check: the gene numbers confirm the FIXED database was used

An earlier Col-CEN `gene.bed` mislabelled between-gene gaps as `intron`
(42.0% intron / 21.6% intergenic / 36.4% exon); it was regenerated on
2026-08-04 and the corrected build gives 20.6 / 43.0 / 36.4.

Col-CEN self-annotation here measures **intron 20.511 / intergenic 42.918 /
exon 36.318** — all three match the corrected figures and none matches the
broken ones. So A3 ran against the fixed `HKS_arabidopsis_ColCEN` v1.2.0, and
the gene-set numbers can be quoted without caveat.
