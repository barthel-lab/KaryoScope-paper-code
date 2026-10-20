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

## Extended Data Fig. 4

`figure/make_build_panels_ab.py` writes panels a and b to
`figures/fig7/build_workflow_panelsAB.svg`; `figure/make_build_figure.py` adds
panel c and writes `figures/fig7/build_workflow.svg`. Panel c reads the two
genome-view karyotype renders in `karyotype_ey15_2/` and `karyotype_colcen/`
(KaryoScope `karyotype` output for Ey15-2 and Col-CEN v1.2 against
`HKS_arabidopsis_ColCEN`, region feature set) and `results/label_composition.tsv`.
Run from this directory:

```bash
python3 figure/make_build_panels_ab.py
python3 figure/make_build_figure.py --a3-dir .
```

Checked 2026-10-20: both outputs are byte-identical to `build_workflow_panelsAB.svg`
and `build_workflow.svg` in KaryoScope-paper PR #35 (`bf6b6c3`), the submitted
Extended Data Fig. 4.
