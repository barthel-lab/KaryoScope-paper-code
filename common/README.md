# common

Code shared by the figure scripts of several analyses.

- `display_names.py`: the display names of the KaryoScope feature labels, the one table every figure
  script uses when it is given `--display-names` (trial, 2026-11-25). A display name is only what a
  figure prints; the labels in the databases and in every data file are unchanged.
  `display_name(label, vocabulary)` takes a label of the v2 database (`KS_human_CHM13_v2`: `active_hor`,
  `hor`, `mon`) or the v1 database (`KS_human_CHM13`: `active_specific`, `inactive_specific`,
  `hor_multigroup1`) and returns, for example, "active HOR", "inactive HOR", "aSat HOR". A parent
  (`<name>_multigroup1`) prints its name alone. `MULTIGROUP` is the text a tree may print after a
  parent row's k-mer count.
- `karyoscope_feature_name_mapping.tsv`: the v1 label of every feature and its v2 label. Every pair
  gets one display name in `display_names.py`, except v1 `noncentromeric_specific` (v2 `rDNA`), which
  prints "noncentromeric".

| figure | script | option |
|---|---|---|
| Fig. 1b (tree) | `hg002/scripts/draw_fig1_trees.py` | `--display-names` (with `--no-multigroup`, `--parent-bars thin`) |
| Fig. 1c (colour key) | `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py` | `--display-names` (before `--key-merge-colours`, `--key-label`) |
| Fig. 3 | `vendored/KaryoScope-ISCN/KaryoScope_assembly_contig_zoom_plot.py` | `--display-names` (`stacked_args: {display_names: true}` in the configuration) |
| Fig. 4, ED Fig. 8, Supp. Fig. 2 | `d4z4/make_figures.py` | `--display-names` |
| ED Fig. 4 | `revision/a3-cross-genome/figure/make_build_figure.py` | `--display-names` |
| ED Fig. 7 | `robs/scripts/figure5_assemble.py` | `--display-names` |

Without the option each script draws what it drew before, byte for byte.
