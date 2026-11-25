# A3 — cross-genome

Reviewer 3 (major point 1) asked whether a database can be built for genomes
other than human. We annotated two non-human genomes with databases built for
them and measured the label composition per feature set.

The per-set breakdown is in `results/label_composition.tsv`, and its columns are
described in `label_composition_README.md` beside it.

The novel fraction is near zero on the genome a database was built from
(Col-CEN nuclear, 0.000% on all four sets) and small but non-zero on a different
accession of the same species (Ey15-2: gene 3.544%, repeat 1.834%, region
0.558%, and chromosome 0.000%). The chromosome set tiles the genome by
construction, and the gene and repeat sets do not.

## Extended Data Fig. 4

`figure/make_build_panels_ab.py` writes panels a and b to
`figures/fig7/build_workflow_panelsAB.svg`. `figure/make_build_figure.py` adds
panel c and writes `figures/fig7/build_workflow.svg`. Panel c reads the two
genome-view karyotype renders in `karyotype_ey15_2/` and `karyotype_colcen/`
(KaryoScope `karyotype` output for Ey15-2 and Col-CEN v1.2 with
`HKS_arabidopsis_ColCEN`, region feature set) and `results/label_composition.tsv`.
Run from this directory:

```bash
mkdir -p figures/fig7
python3 figure/make_build_panels_ab.py
python3 figure/make_build_figure.py --a3-dir . --display-names
```

`--display-names` (2026-11-25) prints the feature labels with the display names of
`../../common/display_names.py` ("p arm", "q arm", "Transposable Element", "centromere gap"; the keys
space their entries by the printed label); without it the labels print as in the database (`p_arm`,
`q_arm`, `Transposable_Element`, `cen_gap`), as checked below.

Checked 2026-10-20: both outputs were identical to `build_workflow_panelsAB.svg`
and `build_workflow.svg` in KaryoScope-paper `main` (`5c11d87`), the submitted
Extended Data Fig. 4. Since 2026-11-21 both scripts set every text in Arial (the
SVGs named Helvetica, then Arial, and a monospace face for the command names
`prep-bed`, `build`, and `annotate` and the priority file in panel b). The
priority file's three columns are now three text elements at the x positions
the monospace setting gave them. Nothing else changed: the same text, shapes,
and positions.

The submitted TIFF is the SVG at 183 mm and 300 ppi, flattened onto white:

```bash
rsvg-convert -w 2161 -b white figures/fig7/build_workflow.svg -o ED_Fig4.png   # then RGB LZW TIFF
rsvg-convert -f pdf -w 183mm figures/fig7/build_workflow.svg -o ED_Fig4.pdf
```
