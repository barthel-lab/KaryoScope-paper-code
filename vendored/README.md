# Vendored code

Only the files that the Fig. 1b-e and Fig. 3 code paths import are copied. The
set was found by running both figures with an import tracer (every module file
loaded by the conductor process and by the figure script it launches). Each
Python file starts with a header naming its source repository, commit and
authors.

## KaryoScope-ISCN (`barthel-lab/KaryoScope-ISCN` at `68a9816`, PR #87)

| file | used by |
|---|---|
| `KaryoScope_assembly_contig_zoom_plot.py` | Fig. 3 (stacked zoom plot, `--panel-letters --hide-scale-bars --no-refine --zoom-levels 5000000,500000,100000`) |
| `KaryoScope_assembly_figure.py` | Fig. 1b-e |
| `KaryoScope_assembly_genome_plot.py` | imported by `KaryoScope_assembly_figure.py` |
| `plot_featureset_tree_bars.py` | imported by `KaryoScope_assembly_figure.py` (hierarchy trees) |

`68a9816` is the head of KaryoScope-ISCN PR #87 (`fig3/lowercase-panel-letters`,
not merged), branched from `main` (`4eb7128`). Its only change to these files is
in `KaryoScope_assembly_contig_zoom_plot.py`: `--panel-letters` draws lowercase
a, b, c (it drew A, B, C), as in the revised manuscript. The other three files
are unchanged from `4e42a98`, the commit the first Fig. 3 render used
(2026-08-28), to `68a9816`. KaryoScope-ISCN is
GPL-3.0-or-later; its `LICENSE` is included.

One change: `plot_featureset_tree_bars.py` reads its resource directory from
`$KS_RESOURCES` when set (the original fixes it to a directory under `$HOME`).
Set it to `vendored/KaryoScope-archive-resources/resources`.

## karyoscope_conductor (`barthel-lab/KaryoScope-conductor` at `0732e96`)

The YAML runner (`python -m karyoscope_conductor --config <yaml>`), at the
commit that KaryoScope-HG002 and KaryoScope-ROBs pin. Copied: the package
`__init__`, `__main__`, `common`, `config`, `pipeline`, `tsv_generator`, and the
two steps the figures use, `steps/iscn/assembly_figure.py` and
`steps/iscn/assembly_contig_zoom_plot_stacked.py`. Steps are imported by name
from the config's `pipeline:` list, so a config that names any other step
fails here. Unmodified.

## KaryoScope-archive-resources

Colour palettes, per-feature k-mer counts (`*.database.info.txt`), priority
files and the v2 hierarchy for the `KS_human_CHM13` (v1) and
`KS_human_CHM13_v2` databases, from the KaryoScope-archive pipeline at commit
`7c976e04` (`resources/`), the commit the HG002 and ROB setups pinned. The later
history of that repository carries an MIT licence (`LICENSE.MIT`).

The figure scripts also import `karyoplot`
(`barthel-lab/KaryoScope-plotlib`, public, GPL-3.0), installed by
`environment/figures.yml` at commit `f59c5fd`.

## Check

On the same inputs, the vendored copy and the original repositories at the
commits above give byte-identical SVG and PNG output, for Fig. 3 and for the
HG002 figure.

## KaryoScope-archive `bin_features.py` (`barthel-lab/KaryoScope-archive` at `4434885`)

| file | used by |
|---|---|
| `KaryoScope-archive/bin_features.py` | binning of the pangenome annotations into the 100 bp (D4Z4) and 100 kb (centromere) windows that `d4z4/` and `centromere/` read |

The script assigns each bin the feature with the largest overlap, ties going to
the first feature seen. It is unchanged below its header. KaryoScope v1.1.0
`karyoscope bin` uses a different rule (a bin is novel only above 50 % novel;
ties broken alphabetically) and does not reproduce these windows.
`revision/a2-k-sensitivity/step4/scripts/v1-region-rebuild/bin_features_feb2026.py`
is the same code under a different header.

