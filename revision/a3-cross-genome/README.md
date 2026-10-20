# A3 — cross-genome

R3-major-1 asks whether a database can be built for genomes other than human.
Two non-human genomes were annotated against purpose-built databases and the
label composition measured per feature set.

`results/label_composition.tsv` carries the per-set breakdown, with
`label_composition_README.md` beside it describing the columns.

Headline: novel fraction is near zero on the genome a database was built from
(Col-CEN nuclear, 0.000% on all four sets) and small but non-zero on a different
accession of the same species (Ey15-2: gene 3.544%, repeat 1.834%, region
0.558%, chromosome 0.000%) — which is the expected direction, since the
chromosome set tiles by construction while the gene and repeat sets do not.
