# HG002 assembly figure (Fig. 1b-e)

`conf/` and `scripts/check_data.sh` are from `barthel-lab/KaryoScope-HG002` at
commit `27885ba`. `requirements.txt` is the original pin list and is kept for
the record; use `../environment/figures.yml` instead. Changes from that commit:

- `conf/hg002.yaml` runs only the `iscn.assembly_figure` step, points
  `scripts.iscn` at `../vendored/KaryoScope-ISCN`, and has a placeholder for
  `colors_dir`;
- `scripts/setup.sh` and `scripts/sync_data.sh` are not copied. They cloned
  private repositories and copied data from the TGen cluster over SSH.

The figure is drawn by `vendored/KaryoScope-ISCN/KaryoScope_assembly_figure.py`
(karyotype grid, plus chromosome, region and subtelomeric hierarchy trees),
called through the vendored `karyoscope_conductor`.

## Inputs

`data/KS_human_CHM13_v2/` must hold the KaryoScope annotation of the HG002
v1.1 diploid assembly with the `KS_human_CHM13_v2` database:

- `HG002.assembly.KS_human_CHM13_v2.<featureset>.presmoothed.scaffolded.binned{100,100000,1000000}.features.bed.gz`
  for `chromosome region repeat gene acrocentric subtelomeric`;
- the centromere BED named in `conf/hg002.yaml`
  (`...centromeres_v1coords.bed`: v1 centromere coordinates mapped onto the v2
  scaffold names).

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

Output: `results/assembly_figure/HG002_figure_v2_v1cen_grid.{svg,png}`.

## Status

The published Fig. 1 is an Illustrator composition built on this render. The
`centromeres_v1coords.bed` file used for it was not found on the cluster, and
the render committed in KaryoScope-HG002 was made with KaryoScope-ISCN
`4d1b086`. Run with the standard
`HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres.bed.gz`
instead, the vendored code produces output byte-identical to the
un-vendored KaryoScope-ISCN `4e42a98` and conductor `0732e96` (the Fig. 1 files are
unchanged from `4e42a98` to the vendored `68a9816`).
