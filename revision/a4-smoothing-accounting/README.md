# A4 — smoothing accounting

What smoothing does, and to how much. Two populations are accounted separately
because they are different problems:

* **novel** — *k*-mers absent from the database, initially labelled with the root
* **multilabel** — *k*-mers whose label is an internal hierarchy node (the LCA of
  several leaves)

## Novel

Before smoothing, **75,396,753 of 2,940,260,101 evaluated positions (2.5643%)**
carry `novel` on HG002 haplotype 1, and this figure is **identical across all six
feature sets**, because novelty is a property of the *k*-mer rather than of the
feature set. That identity is a useful check that the accounting is being read
right.

Run lengths: 99.4% of novel bases lie in runs <= 1 kb, 86.9% in runs <= 31 bp,
and the **median novel run is exactly 31 bp = *k* in every feature set** — a
single substitution destroys up to *k* consecutive *k*-mers, so most novel runs
are point divergence from the reference rather than sequence absent from it. The
455,961 bp in runs longer than 1 kb (0.60% of novel bases, 0.0155% of evaluated
positions) are the candidates for genuinely absent sequence.

## Multilabel

**Multilabel is the dominant ambiguity population and had never been measured**:
pre-smoothing share is 16.03% for the chromosome set against novel's 2.5%.
Smoothing takes chromosome 16.03% → 1.72%.

Per set: chromosome 16.03%, repeat 9.59%, region 8.74%, acrocentric 0.08%,
subtelomeric 0.013%, **gene exactly 0%** — correct, since the gene hierarchy is a
flat three-leaf tiling with no internal nodes. Another good check.

Median resolved run is 3–8 bp everywhere.

## The acrocentric caveat, which is real

Fraction of multilabel bases in runs <= 1 kb, HG002 hap1 / CHM13: subtelomeric
100/100, region 99.5/99.3, repeat 98.8/98.6, chromosome 95.6/93.1, acrocentric
**93.2/67.3**.

**The split is an N-masking artefact, not a CHM13 quirk.** HG002 hap1's N runs
are exactly the acrocentric short arms (chr14 3.46–4.97 Mb, chr15 2.97–5.72,
chr21 2.68–4.98, chr22 2.63–3.71; 7.63 Mb total), and the method masks N runs out
as absence rather than ambiguity — which deletes precisely the rDNA and p-arm
sequence that generates long multilabel runs. CHM13 is gapless so nothing is
masked, and all 427 of its long (>= 1 kb) acrocentric runs are on acrocentrics.

**So 67.3% is the honest number and 93.2% is the flattered one.** This is a real
limitation of `max_gap` in exactly the regions the paper is about, and it is
stated as one.
