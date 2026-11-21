# HG002 assembly figure (Fig. 1b)

`conf/` and `scripts/check_data.sh` are from `barthel-lab/KaryoScope-HG002` at
commit `27885ba`. `requirements.txt` is the original pin list and is kept for
the record. Use `../environment/figures.yml` instead. Changes from that commit:

- `conf/hg002.yaml` runs only the `iscn.assembly_figure` step, points
  `scripts.iscn` at `../vendored/KaryoScope-ISCN`, and has a placeholder for
  `colors_dir`. It has a second analysis, the Fig. 1b render (see Fig. 1).
- `scripts/setup.sh` and `scripts/sync_data.sh` are not copied. They cloned
  private repositories and copied data from the TGen cluster over SSH.

The figure is drawn by `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py`,
called through the vendored `karyoscope_conductor`. The first analysis in
`conf/hg002.yaml` draws the karyotype in three rows of eight chromosomes, plus the
chromosome, region, and subtelomeric hierarchy trees (the render behind the first
submission's Fig. 1b-e). The second draws Fig. 1b.

## Inputs

`data/KS_human_CHM13_v2/` must hold the KaryoScope annotation of the HG002
v1.1 diploid assembly with the `KS_human_CHM13_v2` database:

- `HG002.assembly.KS_human_CHM13_v2.<featureset>.presmoothed.scaffolded.binned{100,100000,1000000}.features.bed.gz`
  for `chromosome region repeat gene acrocentric subtelomeric`.
- the centromere BED named in `conf/hg002.yaml`
  (`...centromeres_v1coords.bed`: v1 centromere coordinates mapped onto the v2
  scaffold names). It is committed here, in
  `data/KS_human_CHM13_v2/` (see Centromere windows).

`scripts/check_data.sh` lists what is missing. On the TGen cluster the BEDs were
under `$KS_PROJECTS/KaryoScope/results/HG002/assembly/combined/KaryoScope/`.

## Run

```bash
conda env create -f ../environment/figures.yml && conda activate ks-paper-figures
export PYTHONNOUSERSITE=1
export KS_RESOURCES=$PWD/../vendored/KaryoScope-archive-resources/resources
sed -i "s|<KS_RESOURCES>|$KS_RESOURCES|" conf/hg002.yaml
PYTHONPATH=../vendored/KaryoScope-conductor python -m karyoscope_conductor --config conf/hg002.yaml
```

Output: `results/assembly_figure/HG002_figure_v2_v1cen_grid.{svg,png}` and
`results/assembly_figure/HG002_fig1b_two_rows.{svg,png}`.

## Fig. 1

Fig. 1 has two panels at the Nature main-figure size (180 mm wide, at most 170 mm
tall): a, the FISH and KaryoScope schematic, and b, the HG002 karyotype. The
hierarchy trees of the first submission's Fig. 1c-e are Extended Data Figs. 1b,
2b, and 3b.

Panel b is the second analysis in `conf/hg002.yaml`, which runs
`KaryoScope_assembly_figure.py` with these options, added in this repository
(without them the output is unchanged):

| option | effect |
|---|---|
| `--chroms-per-row 12` | two rows: chr1 to chr12, then chr13 to chr22, X, and Y |
| `--no-trees` | no hierarchy trees (the scale legend goes below the last row) |
| `--no-legend-colons` | legend rows read Genome, Centromere, and Telomere, as in the published panel |
| `--width-mm 180` | page written 180 mm wide, with text sized in printed points: chromosome names 6.5 pt (`--chrom-pt`), h1 and h2 5.5 pt (`--hap-pt`), Mb labels 5 pt (`--size-pt`), and legend 6 pt (`--legend-pt`) |
| `--margin 6` | 6-unit page margin |

The zoom factors do not depend on the arrangement: the legend prints 4.9x for the
centromere windows and 2158x for the 25 kb telomere windows, as in the first
submission. The inputs and the centromere window file are the same.

`scripts/compose_fig1.py` places panel a above panel b and writes the figure. Panel
a is cropped from the published Illustrator file
(`KaryoScope-paper/figures/fig1/HG002_figure_v2_v1cen_grid.ai`, PDF-compatible, not
modified) as vector artwork. Everything below panel a is deleted from the copy
first. Panel a is scaled so its smallest text prints at 5.05 pt (141 x 57 mm), and
its row labels FISH and KaryoScope are set again in Arial Bold at 7 pt. Panel
letters are 8 pt bold lowercase Arial.

The script needs PyMuPDF and cairosvg (`environment/d4z4.yml` has cairosvg; add
PyMuPDF with `pip install pymupdf==1.26.5`).

```bash
# fontconfig must return Arial for "sans-serif" (cairosvg converts the SVG text)
python scripts/compose_fig1.py \
    --ai /path/to/KaryoScope-paper/figures/fig1/HG002_figure_v2_v1cen_grid.ai \
    --panel-b results/assembly_figure/HG002_fig1b_two_rows.svg \
    --arial-bold /path/to/arialbd.ttf --out results/Fig1
```

Output: `results/Fig1.pdf` (vector, fonts embedded) and `results/Fig1.png` (300
ppi), 180 x 164.5 mm, text 5.0 to 7.0 pt. The script stops if the PDF holds a
font other than Arial or if the page is taller than 170 mm.

## Centromere windows (Fig. 1b)

The Fig. 1b zoom windows are the rows of
`data/KS_human_CHM13_v2/HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres_v1coords.bed`
(46 rows, sha256 `b7390b9381d49a6c3f0f4e4e4583f8813aefc2c00c1fe97da1af7b2a33d29080`).
The original file was not found. This one was rebuilt with
`scripts/make_centromeres_v1coords.py` from the two scaffolded centromere BEDs
of the same assembly:

| input | database | sha256 |
|---|---|---|
| `HG002.assembly.KS_human_CHM13.region.presmoothed.scaffolded.centromeres.bed` | `KS_human_CHM13` (v1) | `151b0a121cbfd92c08fa0ee8397707828dcef7a162ae41856c0963e798502754` |
| `HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres.bed.gz` | `KS_human_CHM13_v2` | `8243e6bd036cd27474f56cc9a2de878458e3f58ecefe634d9c3ac6855c6eac64` (uncompressed `0e9a7c26def1d979024bbf0a4e391874c7f0fd88d19a8fbb625b8f0cdf18fb63`) |

```bash
python3 scripts/make_centromeres_v1coords.py <v1 centromeres.bed> <v2 centromeres.bed.gz> \
    > data/KS_human_CHM13_v2/HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres_v1coords.bed
```

The two inputs list the same 46 chromosomes in the same order. Each output row
takes the v2 scaffold name (`chr1_hap1_chr1_PATERNAL`) and the v1 coordinates
(row of `chr1_hap1_chr1_hap1`). All 46 window lengths equal the 46 zoom-window
lengths printed in Fig. 1b.

## Status

The first submission's Fig. 1 is an Illustrator composition built on the
three-row render. It was exported with `illustrator/nature_ai_export.py --relabel C=E,E=C`, so panel c is the chromosome tree and panel e the subtelomeric tree. The
render committed in KaryoScope-HG002 was made with KaryoScope-ISCN `4d1b086`.
Run with the standard
`HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres.bed.gz`
instead of the window file, the vendored code produces output identical
to that of the un-vendored KaryoScope-ISCN `4e42a98` and conductor `0732e96` (the
Fig. 1 files are unchanged from `4e42a98` to the vendored `68a9816`).
