# A4 — smoothing accounting

We measured what smoothing does, and to how much sequence. We accounted for two
populations separately:

* **novel**: *k*-mers absent from the database, initially labelled with the root
* **multilabel**: *k*-mers whose label is an internal hierarchy node (the lowest
  common ancestor, LCA, of several leaves)

## Novel

Before smoothing, **75,396,753 of 2,940,260,101 evaluated positions (2.5643%)**
on HG002 haplotype 1 are labelled `novel`, and this figure is **identical across
all six feature sets** because they were all built from the same set of
*k*-mers.

Run lengths: 99.4% of novel bases lie in runs <= 1 kb, 86.9% in runs <= 31 bp,
and the **median novel run is exactly 31 bp = *k* in every feature set**. A
single substitution changes up to *k* consecutive *k*-mers, so most novel runs
are point differences from the reference rather than sequence absent from it.
The 455,961 bp in runs longer than 1 kb (0.60% of novel bases, 0.0155% of
evaluated positions) are the candidates for absent sequence.

## Multilabel

**Multilabel positions are the larger ambiguous population, and had not been
measured before.** Before smoothing they are 16.03% of the chromosome set,
compared with 2.5% novel. Smoothing reduces the chromosome-set share from
16.03% to 1.72%.

Per set: chromosome 16.03%, repeat 9.59%, region 8.74%, acrocentric 0.08%,
subtelomeric 0.013%, and **gene exactly 0%**. The gene hierarchy is a flat tiling
of three leaves with no internal nodes, so 0% is expected.

The median resolved run is 3–8 bp in every feature set.

## The acrocentric caveat

Fraction of multilabel bases in runs <= 1 kb, HG002 hap1 / CHM13: subtelomeric
100/100, region 99.5/99.3, repeat 98.8/98.6, chromosome 95.6/93.1, acrocentric
**93.2/67.3**.

**The difference comes from N masking, not from CHM13.** The N runs of HG002
hap1 are exactly the acrocentric short arms (chr14 3.46–4.97 Mb, chr15 2.97–5.72,
chr21 2.68–4.98, and chr22 2.63–3.71, 7.63 Mb in total). N runs are masked as
absent sequence rather than ambiguity, which removes the rDNA and p-arm sequence
that produces long multilabel runs. CHM13 has no gaps, so nothing is masked, and
all 427 of its long (>= 1 kb) acrocentric runs are on acrocentric chromosomes.

**The CHM13 value (67.3%) is the one to quote. The HG002 value (93.2%) is
inflated by the masking.** This is a limitation of `max_gap` in the regions the
paper studies, and we state it as one.
