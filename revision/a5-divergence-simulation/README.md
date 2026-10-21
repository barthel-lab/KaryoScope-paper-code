# A5 — divergence simulation

We measured how far a query can diverge from the reference before annotation
degrades. `mutate_fasta.py` introduces substitutions at a set rate. All 13
conditions ran to completion with 100% coverage.

The driver is the concordance harness in `../concordance`, so the numbers are
comparable with the accuracy figures in the manuscript.
