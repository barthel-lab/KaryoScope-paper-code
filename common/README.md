# common

Code shared by the figure scripts of several analyses.

- `display_names.py`: the display names of the KaryoScope feature labels, the one table every figure
  script uses when it is given `--display-names` (2026-11-25). A display name is only what a
  figure prints; the labels in the databases and in every data file are unchanged.
  `display_name(label, vocabulary)` takes a label of the v2 database (`KS_human_CHM13_v2`: `active_hor`,
  `hor`, `mon`) or the v1 database (`KS_human_CHM13`: `active_specific`, `inactive_specific`,
  `hor_multigroup1`) and returns, for example, "active HOR", "inactive HOR", "aSat HOR". A parent
  (`<name>_multigroup1`) prints its name alone. `MULTIGROUP` is the text a tree may print after a
  parent row's k-mer count. `display_sort_key(name, featureset)` orders legends in hierarchy
  order: the rows of the feature set's tree from top to bottom (the v2 trees of Fig. 1b and ED Figs. 2-3,
  children by priority), so a parent is listed directly above its children -- for region centromeric,
  aSat, aSat HOR, active HOR, divergent HOR, inactive HOR, mixed aSat, monomeric, bSat, cenSat, ct, gSat,
  HSat, HSat1, HSat1A, HSat1B, HSat2, HSat3, rDNA, arm, p arm, q arm -- then any other name
  alphabetically. Each legend section passes its feature set (rDNA is a row of both the region and the
  acrocentric tree); the catch-all and "novel" entries stay last. `display_text` rewrites a printed phrase (its feature
  labels as display names, beta-satellite as bSat). `shared_parent_name` names a swatch several parents
  share: the deepest of them in the v1 hierarchy, the root excluded ("acrocentric").
- `karyoscope_feature_name_mapping.tsv`: the v1 label of every feature and its v2 label. Every pair
  gets one display name in `display_names.py` (v1 `noncentromeric_specific` is v2 `rDNA` and prints
  "rDNA").

| figure | script | option |
|---|---|---|
| Fig. 1b (tree) | `hg002/scripts/draw_fig1_trees.py` | `--display-names --no-multigroup` (`--parent-bars thin` is an alternative) |
| Fig. 1c (colour key) | `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py` | `display_names: true` in `hg002/conf/hg002.yaml` (applied before `--key-merge-colours`, `--key-label`) |
| Fig. 3 | `vendored/KaryoScope-ISCN/KaryoScope_assembly_contig_zoom_plot.py` | `--display-names` (`stacked_args: {display_names: true}` in the configuration) |
| Fig. 4, ED Fig. 8, Supp. Fig. 2 | `d4z4/make_figures.py` | `--display-names` (legends, "Other Elements" key, ED Fig. 8a schematic) |
| ED Fig. 4 | `revision/a3-cross-genome/figure/make_build_figure.py` | `--display-names` |
| ED Fig. 7 | `robs/scripts/figure5_assemble.py` | `--display-names` |

The documented commands (READMEs, `hg002/conf/hg002.yaml`, `robs/conf/config.template.yaml`) use the
option. Without it each script draws the submitted figure, byte for byte.
