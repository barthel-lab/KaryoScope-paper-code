# Robertsonian translocations (Fig. 3, Extended Data Fig. 7)

Code from `barthel-lab/KaryoScope-ROBs` at `1188850`, the head of PR #2
(`nature/ed-spec`, on top of branch `fig5/assemble-panels`; not merged), which
sets Extended Data Fig. 7 at the Nature size (183 mm wide, 5-7 pt text, 8 pt
lowercase panel letters). Only the files the two figures use are copied:

| file | role |
|---|---|
| `conf/config.template.yaml` | the figure configuration, with every per-sample value replaced by a placeholder; its `a2_step3`, `a6_rob` and `a7_job3` sections are read by the revision analyses in `../revision/` |
| `scripts/breakpoint_density_stacked.py` | Extended Data Fig. 7 a-c (breakpoint refinement per sample) |
| `scripts/sst1_kmer_composition.py` | Extended Data Fig. 7 d (chromosome-specific k-mers in the SST1 arrays) |
| `scripts/chm13_kmer_lollipop.py` | Extended Data Fig. 7 e (positions of those k-mers in T2T-CHM13; CHM13 coordinates) |
| `scripts/figure5_assemble.py` | composes Extended Data Fig. 7 from the panels |
| `scripts/annotate_scaffold_v110.py` | step 2 a: KaryoScope v1.1.0 annotate and scaffold, presmoothed (new here) |
| `scripts/check_data.sh` | lists missing input BEDs |
| `requirements.txt` | the original pin list; use `environment/figures.yml` |

Fig. 3 is drawn by `vendored/KaryoScope-ISCN/KaryoScope_assembly_contig_zoom_plot.py`
through the vendored conductor (`vendored/KaryoScope-conductor`), from the
`assembly_contig_zoom_plot_stacked_analyses` section of the configuration.

Changes from `1188850`:

- `breakpoint_density_stacked.py` reads the configuration from `$ROBS_CONFIG`
  when set (default `conf/robs.yaml`), and its `--db-root` default is the
  registry under `$KS_SOFTWARE` instead of a cluster path;
- `sst1_kmer_composition.py` reads its SST1 array table from the configuration's
  `sst1_arrays` section; at `1188850` the table is written in the script;
- `check_data.sh`: comments and the failure message point at `data/robs/`
  instead of the cluster sync script, which is not copied;
- comments in `breakpoint_density_stacked.py` and `figure5_assemble.py` no
  longer quote per-sample bin counts or display-window widths;
- `scripts/annotate_scaffold_v110.py` is new (step 2 a);
- `conf/robs.yaml` is replaced by the template, and `scripts/FIGURE5_REPRODUCE.md`
  is not copied (both carry per-sample coordinates). The steps are below.

## Steps

The commands use placeholders: `<asm>` is the assembly name used in the file
names (`GM03417.assembly`, `GM03786.assembly`, `GM04890.assembly`), `<scaffold>`
a fusion scaffold, `<root>` the KaryoScope database root.

### 1. Data access

The three cell-line assemblies (GM03417, GM03786, GM04890; two haplotype FASTAs
each) are controlled-access: request access to dbGaP study **phs003920.v1.p1**.
Keep the assemblies, every annotation made from them and the filled
configuration under the terms of that access.

### 2. Annotate, scaffold and bin

All three steps were checked on one sample (see What reproduces). Versions:
**annotate and scaffold with KaryoScope v1.1.0** (tag `4c09ce5`,
doi:10.5281/zenodo.20657817) and the `KS_human_CHM13` (v1, KMC) database from
the KaryoScope database registry; **bin with `karyoscope bin`** from v1.1.0 or
v2.3.0 (both give the same figures).

a. Annotate and scaffold both haplotypes of each assembly (about 70 min and
40 GB per assembly on 16 threads; KMC needs the memory):

```bash
# KaryoScope v1.1.0 on PYTHONPATH (or installed), its get_featureIDs built (native/get_featureIDs: make)
export KARYOSCOPE_GET_FEATUREIDS=<KaryoScope v1.1.0>/native/get_featureIDs/build/get_featureIDs
python scripts/annotate_scaffold_v110.py --hap1 <hap1 FASTA> --hap2 <hap2 FASTA> \
    --db-root <root> -o <annotation_dir> -t 16
```

The script calls v1.1.0's scaffold pipeline in BED mode on the presmoothed
annotations; v1.1.0's `karyoscope scaffold` command line writes only smoothed
scaffolded BEDs, and `karyoscope karyotype --presmoothed` stops on the v1
database's colour table. Scaffold names are `<chrom>_<hap>_<contig>[_rc]`, with
the v1 chromosome label (for example `chr14_specific_hap1_<contig>`).

b. Join the two haplotypes into the per-assembly per-base BEDs the figures read
(KaryoScope output is in genome order; do not sort):

```bash
mkdir -p data/robs
for set in chromosome acrocentric region repeat gene subtelomeric; do
  cat <annotation_dir>/<hap1 stem>.KS_human_CHM13.$set.presmoothed.scaffolded.bed.gz \
      <annotation_dir>/<hap2 stem>.KS_human_CHM13.$set.presmoothed.scaffolded.bed.gz \
      > data/robs/<asm>.KS_human_CHM13.$set.presmoothed.scaffolded.bed.gz
done
```

c. Bin at 100 bp, 100 kb and 1 Mb, without `--db` (plain largest-overlap
bins). Only these three sizes: the Fig. 3 zoom plot uses the largest available
bin that divides each window, so extra sizes change the figure. Bins are per
scaffold, so binning only the fusion scaffold's rows gives the same bins. The
output name must end in `.bed.gz` (`karyoscope bin` picks the format from the
extension):

```bash
for set in chromosome acrocentric region repeat gene subtelomeric; do
  for b in 100 100000 1000000; do
    karyoscope bin -i data/robs/<asm>.KS_human_CHM13.$set.presmoothed.scaffolded.bed.gz -b $b -t 8 \
        -o data/robs/<asm>.KS_human_CHM13.$set.presmoothed.scaffolded.binned$b.features.bed.gz
  done
done
```

The KaryoScope-archive `bin_features.py` (as used for the D4Z4 and centromere
windows) does not give these figures: on the test sample, 41 of the 1,056 drawn
rectangles of its Fig. 3 panel change (30 drawn in their place), mostly
chromosome multigroup.

d. For Extended Data Fig. 7 d, a BGZF copy of the chromosome BED with a tabix
index:

```bash
zcat data/robs/<asm>.KS_human_CHM13.chromosome.presmoothed.scaffolded.bed.gz | bgzip -c \
    > data/robs/<asm>.KS_human_CHM13.chromosome.presmoothed.scaffolded.features.bed.gz
tabix -p bed data/robs/<asm>.KS_human_CHM13.chromosome.presmoothed.scaffolded.features.bed.gz
```

`scripts/check_data.sh` lists anything still missing.

### 3. Find the fusion scaffold and the SST1 block midpoint

Methods: the Fig. 3 zoom windows are the full contig and windows of 10 Mb,
1 Mb and 200 kb centred on the midpoint of the SST1 block at the fusion; the
Extended Data Fig. 7 display window is chosen by hand around each fusion and
extended by 1 Mb on each side.

- **Fusion scaffold.** In each assembly one scaffold joins the long arms of the
  two acrocentric chromosomes of the translocation (rob(14;21) in GM03417,
  rob(13;14) in GM03786 and GM04890). In the chromosome annotation its labels
  change from one of the two chromosomes to the other along the scaffold, and
  in the acrocentric annotation it carries a run of `SST1_specific` between the
  two q arms. The KaryoScope-ISCN report lists it as a breakpoint. Take its name
  from the scaffold map (`scaffold`) and the contig name from the assembly
  (`original_scaffold`).
- **SST1 block midpoint.** Take the `SST1_specific` rows of the fusion scaffold
  from the per-base acrocentric BED and merge runs closer than 1 kb:

  ```bash
  zcat data/robs/<asm>.KS_human_CHM13.acrocentric.presmoothed.scaffolded.bed.gz \
    | awk -F'\t' -v s=<scaffold> '$1==s && $4=="SST1_specific"' \
    | bedtools merge -d 1000 -i - | awk -v OFS='\t' '{print $0, $3-$2}'
  ```

  The longest merged run is the main SST1 array (the same run for merge
  distances from 500 bp to 10 kb). Round its start and end to the nearest
  100 bp; their midpoint is `center` of the `*_robertsonian` group. Fig. 3 is
  drawn with `--no-refine`, so the zoom windows are centred exactly there.
- **Partial arrays (Extended Data Fig. 7 d).** The upstream and downstream
  partial arrays are the nearest merged `SST1_specific` runs of at least 500 bp
  on each side of the main array. Their bounds in the original configuration
  were set by hand to about 100 bp.

On the test sample these rules, applied to fresh v1.1.0 annotations, give the
configured fusion scaffold (up to the `_specific` label in its name), the
configured main-array bounds and `center` exactly, and the partial arrays within
about 100 bp.

### 4. Fill the template

Copy `conf/config.template.yaml` to `conf/robs.yaml` (or anywhere, and set
`ROBS_CONFIG` to it) and replace every `<...>`:

- `inputs`: the assembly and annotation paths (a record; no script reads it);
- `scripts.iscn`: the path to `vendored/KaryoScope-ISCN`; `colors_dir`:
  `vendored/KaryoScope-archive-resources/resources/databases/KS_human_CHM13/`;
- `breakpoint_groups.*_robertsonian` (Fig. 3 and Extended Data Fig. 7):
  `scaffold` and `original_scaffold` from step 3; `center` = the SST1 block
  midpoint; `start` and `end` = the display window, chosen by hand around the
  fusion (Extended Data Fig. 7 a-c add 1 Mb on each side; they are not a
  confidence interval);
- `breakpoint_groups.*_all` and `*_breakpoints_only`: the KaryoScope-ISCN report
  rows for each sample (not needed for the two figures);
- `sst1_arrays` (Extended Data Fig. 7 d): the scaffold and, per sample, the
  start and end of the upstream partial array, the main array and the
  downstream partial array from step 3.

### 5. Run

In `environment/figures.yml`, from this directory, with `data/robs/` filled as
in step 2. Extended Data Fig. 7 d also needs `tabix` (htslib) on `PATH`. The
figures were rendered with Arial: install it and let fontconfig resolve
`sans-serif` and `Helvetica` to it, or the text falls back to DejaVu Sans and
the layout shifts.

```bash
conda activate ks-paper-figures; export PYTHONNOUSERSITE=1 ROBS_CONFIG=conf/robs.yaml
C=../vendored/KaryoScope-archive-resources/resources/databases/KS_human_CHM13
# Fig. 3 -> results/contig_zoom_plots/robertsonian_integrated.svg
PYTHONPATH=../vendored/KaryoScope-conductor python -m karyoscope_conductor --config $ROBS_CONFIG
python scripts/fig3_headroom.py results/contig_zoom_plots/robertsonian_integrated.svg out/Fig3.svg   # submitted Fig. 3
# Extended Data Fig. 7 a-c, then the a-e composition
python scripts/breakpoint_density_stacked.py --resolution 100 --colors-dir $C --no-legend \
    --emit-anchors out/panelABC_anchors.json --output-dir out/abc
python scripts/figure5_assemble.py --abc-svg out/abc/breakpoint_density_stacked.svg \
    --abc-anchors out/panelABC_anchors.json --bed-dir data/robs --colors-dir $C \
    --region-bed <CHM13 v1 region presmoothed BED> --chroms chr7,chr9,chr17,chr20 \
    --output out/figure5_assembled.svg
```

`--db-root` (default `$KS_DB_ROOT`, else `~/.karyoscope/db`) is read only for legend grouping, which `--no-legend`
turns off. Panels a-c take about 8 min and 19 GB. The CHM13 region BED for panel e is
T2T-CHM13 v2.0 annotated with the same database:
`karyoscope annotate -i chm13v2.0.fasta -o <dir> --db KS_human_CHM13 --db-root <root> --feature-set region -t 16`
(the `.presmoothed.bed.gz` output).

## What reproduces

Tested 2026-10-20 on one sample (the assembly of Fig. 3 a), from its two
haplotype FASTAs, with the commands above (step 2 a run as
`scripts/annotate_scaffold_v110.py`; bins from v1.1.0 `karyoscope bin`):

- **Annotation:** the v1.1.0 presmoothed BEDs of both haplotypes are identical
  to the ones the published inputs were made from (KaryoScope-archive pipeline),
  in all six feature sets (smoothed: identical once the archive's trailing
  `_specific` strip is applied).
- **Scaffolding:** v1.1.0 places and orients 59 of the 61 contigs as the
  archive pipeline did, the fusion scaffold among them; two contigs of
  haplotype 2, not drawn in either figure, differ. Names differ only by the
  `_specific` in the chromosome label, which is not drawn. The per-base rows of
  the fusion scaffold are identical in all six feature sets.
- **Fig. 3:** rendered through `vendored/KaryoScope-conductor` and
  `vendored/KaryoScope-ISCN` with these inputs for that sample (and, for the other
  two samples, the archive pipeline's per-base BEDs binned the same way), the SVG is
  byte-identical to the submitted Fig. 3 before one edit: all 3,250 drawn shapes
  and 92 texts. The submitted SVG then has 4 units of headroom added (the top scale
  labels were drawn above the canvas): `viewBox="0 -4 647 809"`, height 809, and
  the background rect at y = -4. `scripts/fig3_headroom.py <render.svg> <out.svg>`
  applies it; render plus this edit equals KaryoScope-paper `main` (`5c11d87`)
  `figures/fig4/robertsonian_integrated.svg` byte for byte.
  The published panel's own annotation files were not kept; this reproduces
  its data and its pixels.
- **Extended Data Fig. 7:** these scripts with the filled configuration
  reproduce panels a-c and the a-e composition byte for byte (SVG, anchors and
  `breakpoint_results.tsv`) from the archive pipeline's per-base BEDs, with the
  100 bp bins made by `karyoscope bin`. The test sample's per-base inputs from
  fresh v1.1.0 annotations are identical to those, so its panels are too.

Binning is the step that matters. The KaryoScope-archive pipeline that first
annotated these assemblies also wrote binned files, made with `bin_features.py`
(on the test sample byte-identical to a rerun of the vendored copy); drawn from
those, Fig. 3 has 3,232 shapes, 117 rectangles differing from the submitted
figure. The same pipeline's per-base BEDs binned with `karyoscope bin` (v1.1.0
or v2.3.0) give the submitted figure.
