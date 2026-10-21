# Centromere clustering (Fig. 5a, 5c, 5d)

`scripts/` is `barthel-lab/KaryoScope-analysis`, branch `pangenome_structure_V2`,
commit `90e2874`, directory `scripts/`, with three changes. The default
`KARYOSCOPE_DIR` in `HPRC2_centromere_analysis.sh` is now a placeholder, so set
the variables below. A font path in `KaryoScope_plot_reads.py`, which is not on
the figure path, is now `$KS_FONT_FILE`. A read name used as a docstring example
in `KaryoScope_cluster_plot.py` and `KaryoScope_visualize_translocation_reads.py`
is replaced by a placeholder.

`scripts/HPRC2_centromere_analysis.sh` runs the whole analysis:

| step | script | output |
|---|---|---|
| 1 | `concatenate_centromere_regions.R` | the 24 per-chromosome BEDs as one BED |
| 2 | `QCfilter_apply.R` (`Err=0 COLLAPSE=0 COLLAPSE_VAR=0`) | `*.pass.bed` |
| 3 | `KaryoScope_cluster_analysis.py` (Ward, structure mode, k by silhouette) | `allchr_structure.sequence_assignments.tsv` |
| 4 | `KS_allchr_dendrogram.py --sil-threshold 0.5 --centroid-sd 5` | `allchr_dendrogram.svg` |
| 5 | `KS_allchr_annotate.py` | `allchr_outlier_annotations.tsv` (48 rows) |
| 6 | `KS_allchr_dendrogram.py --annotations ...` | `allchr_dendrogram_annotated.svg` (Fig. 5a) |
| 7a | `KS_allchr_barplot.py` | `allchr_barplot.svg` (Fig. 5c), and prints the Major/Outlier counts |
| 7b | `KS_allchr_allele_heatmap.py` | `allchr_allele_heatmap*.svg` and the allele tables (Fig. 5d) |

Steps 4-7 apply the two-stage gate: chromosomes whose k = 2 silhouette is
below 0.5 are collapsed to one cluster, then major haplotypes more than 5
standard deviations from their cluster centroid become `chrN_Outlier_S2_i`. The acrocentric
chromosomes and chrY are clustered in step 3 and dropped at
`KS_allchr_dendrogram.py` (`EXCLUDED_CHROMS`).

## Run

```bash
conda env create -f ../environment/centromere.yml     # ks-paper-centromere
conda activate ks-paper-centromere
export PYTHONNOUSERSITE=1
PER_CHR_BED_DIR=<dir with pangenome.chr{1..22,X,Y}.centromere.KS_human_CHM13.presmoothed.region.bed> \
QC_TSV=<NucFlag/Flagger QC table, see below> \
COLORS_DIR=<KS_human_CHM13 database directory, for the colour palettes> \
OUT_DIR=<output dir> \
  bash scripts/HPRC2_centromere_analysis.sh
```

About 30 s on 4 cores, with peak memory under 64 GB.

Set `PYTHONNOUSERSITE=1`, because otherwise packages in `~/.local` shadow the
environment. The result depends on the scipy and scikit-learn versions, and the
versions that give the published counts are pinned in the environment file.

## Inputs

- **Per-chromosome centromere BEDs**: KaryoScope `region` feature set, v1
  database `KS_human_CHM13`, pre-smoothing, on the HPRC Release 2 assemblies,
  clipped to each centromere interval and binned at 100 kb. How they were made
  is described in the manuscript Methods. On the TGen cluster they are
  `$KS_PROJECTS/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed/`
  (93,452 rows, 6,138 contigs).
- **QC table**: `qc/tables/aggregate_qc_v4.curated.tsv`, one row per
  (sample, haplotype, chromosome, region, contig) with the per-region counts of
  each Flagger/NucFlag flag (`_n_Col`, `_n_Dup`, `_n_Err`, `_n_NNN`, `_n_COLLAPSE`,
  `_n_COLLAPSE_OTHER`, `_n_COLLAPSE_VAR`, `_n_MISJOIN`, `_n_HET`). 18,414 data
  rows, 6,138 of them `centromere`. It is rebuilt in `qc/` from the HPRC Release 2
  QC calls (HMM-Flagger v1.2.0 conservative calls in their per-haplotype,
  gap-masked form, and NucFlag v0.3.3), where each step is described. The run used
  Yi-An Chen's copy, `aggregate_qc_v4.manual_curation (3).tsv` (sha256
  `e149eed58c72cfecc832f9feaf93455fc6b23c66cadcdc636cd41f6cc033f64e`), which
  differs from the committed table only in line endings (CRLF). Both give the
  same `*.pass.bed`.
- **Colour palettes**: these affect only the SVG colours.

## Expected output

Checked on 2026-10-20 with `environment/centromere.yml` (Python 3.12.13,
numpy 2.4.3, pandas 3.0.1, scipy 1.17.1, scikit-learn 1.8.0, matplotlib 3.10.8,
R 4.4.1):

- `*.pass.bed` sha256 `eefde8493f147e01db01723020a56555c5f6d6318a3ac497ce75386f053b6e06`.
- Step 3 keeps 4,646 of the 4,650 QC-passed centromeres (four have intervals
  longer than `--max-sequence-length 50000000`).
- Step 7a prints `Grand total: 4551 (3730 Major, 821 Outlier)`.
- `allchr_outlier_annotations.tsv` is identical to
  `revision/a6-parameter-sweeps/centromere-clustering/results/allchr_outlier_annotations.tsv`
  (48 rows, sha256 `d5be5c31…097bd0`, with the `n` column summing to 821).
- `allchr_allele_summary.tsv`, `allchr_allele_pairs.tsv`, and
  `allchr_allele_cooccurrence.tsv` are identical to the copies committed
  on `pangenome_structure_V2` at `90e2874` (`agent_results/`).

The gated per-haplotype assignments are held in memory only and are not
written to disk.

## Note: one chromosome passed QC without being evaluated

HG02965 hap1 chr14 (`HG02965#1#CM087353.1`) is chromosome-level, but its three
QC regions are empty (the step that wrote the regions did not treat it as
chromosome-level, while the step that wrote the table did). Its rows in the QC
table are therefore all `.`, which counts as zero flags, and it passes. It is
in the 4,650 QC-passed and the 4,646 clustered centromeres, and not in the
4,551 analysed, because chr14 is excluded downstream. Over its clustering
interval it has 1 `COLLAPSE_VAR` call. The QC filter is left as it ran, so the
4,650 is what the pipeline produced.
