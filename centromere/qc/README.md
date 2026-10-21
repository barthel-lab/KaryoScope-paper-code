# Centromere assembly QC table (input to `QCfilter_apply.R`)

This directory rebuilds the QC table that `../scripts/QCfilter_apply.R` reads
(step 2 of `../scripts/HPRC2_centromere_analysis.sh`), from the public HPRC
Release 2 QC calls. With `Err=0 COLLAPSE=0 COLLAPSE_VAR=0` it keeps 4,650 of the
6,138 chromosome-level centromeres (ED Fig. 9a "Unfiltered 6,138").

## The chain

| step | what | here | output |
|---|---|---|---|
| 1 | count every QC call over three regions of each chromosome-level contig | `build_master_qc_table.py` | `tables/master_qc_table.full.tsv` |
| 2 | count nine flags per (sample, haplotype, chromosome, region, contig), tools merged | `aggregate_qc.py` | `aggregate_qc_v4.tsv` (not committed) |
| 3 | drop the 2 rows of one non-chromosome-length contig (below) | `aggregate_qc.py --curate` | `tables/aggregate_qc_v4.curated.tsv` |
| 4 | keep contigs with `Err=0 COLLAPSE=0 COLLAPSE_VAR=0` | `../scripts/QCfilter_apply.R` | `*.pass.bed` (4,650 contigs) |

Steps 1 and 2 first ran in the KaryoScope-archive Snakemake pipeline
(2025-11-12 to 2025-11-19) and as an awk script; `original/` holds those files
as a record (see Provenance). `build_master_qc_table.py` and `aggregate_qc.py`
are standalone versions that need only Python 3 (standard library). They give
the same bytes.

### Inputs to step 1

- **QC calls**: the HPRC Release 2 assembly QC calls, one BED per haplotype and
  tool (`prerelease_qc_index.tsv`, 1,386 files, 462 per tool):
  - `flagger_hifi`, `flagger_ont`: HMM-Flagger v1.2.0 conservative calls in
    their per-haplotype, gap-masked form,
    `*.hmm_flagger_v1.2.0_prediction.conservative.canonical.no_Hap.hap{1,2}.bed`.
    In these files assembly gaps are labelled `NNN`. They sit on S3 next to the
    released `*.conservative.bed`, which has no gap label and adds `Hap` calls;
    the two give different counts (see `sensitivity/`).
  - `nucflag`: NucFlag v0.3.3, `{sample}_{hap1,hap2,pat,mat}_hprc_r2_v1.x.y.nucflag.bed`.
  - Index files the URLs come from: `flagger_hifi_v0.1.csv`, `flagger_ont_v0.1.1.csv`,
    `nucflag_pre_release_v0.1.index.csv`, and, for the 6 NucFlag files of
    HG01978, HG02257 and HG03516 (v1.1.0 assemblies),
    `nucflag_pre_release_v0.2.index.csv`. Column `index_csv` names the index of
    each file.
  - The BEDs are not in this repository. `download_prerelease_qc.sh <dir>`
    fetches them (public bucket; AWS CLI without credentials, or curl) into
    `<dir>/<tool>/<file>` and checks the sha256 of each. The sha256 values are
    those of the files the 2025 run read.
- **`inputs/cm_contigs.tsv`**: the 6,139 chromosome-level contigs (GenBank `CM`
  accessions) and the chromosome KaryoScope assigned to each, from the v1
  database (`human_assembly_haploid`) annotation at 1 Mb, with contig lengths
  from the Release 2 `.fai`. 6,138 (sample, haplotype, chromosome)
  combinations; one lists two contigs (below).
- **`inputs/regions.tsv`**: the three regions of each chromosome-level contig,
  0-based half-open, 18,411 rows (6,137 contigs x 3):
  - `subtelo_p` = `[0, 500,000)`;
  - `subtelo_q` = `[L - 500,000, L)`, L = contig length;
  - `centromere` = the span, from the first start to the last stop, of the
    1 Mb bins labelled centromeric in an earlier KaryoScope annotation of the
    same assembly (v1 database `human_assembly_haploid`, 1 Mb bins,
    `*.featureIDs.binsize1000000.centromeres.scaffolded.bed`). No contiguity
    test: every centromere-class bin counts.

  **This file is the only surviving record of the centromere regions.** The 1 Mb
  annotation files they were computed from were deleted, so the centromere
  regions cannot be re-derived; `subtelo_p` and `subtelo_q` can (from the
  `.fai`).

### Conventions of the master table

Rows nest as sample (order of `prerelease_qc_index.tsv`, HG005 first) x
haplotype x chromosome (chr1-22, X, Y) x region (`subtelo_p`, `centromere`,
`subtelo_q`) x tool (`flagger_hifi`, `flagger_ont`, `nucflag`):

- not chromosome-level: one row `N/A N/A N/A N/A`;
- chromosome-level, no call in the region: one row `<contig> . . .`;
- otherwise one row per call with at least 1 bp in the region, in QC-file
  order, clipped to the region (as `bedtools intersect` without `-wa`). 770 of
  the 159,002 call rows are clipped.

Step 2 counts `Col Dup Err NNN COLLAPSE COLLAPSE_OTHER COLLAPSE_VAR MISJOIN HET`
per key and keeps the awk's column names (`_n_Col`, ...; R reads them as
`X_n_Col`, ...). `.` counts as zero. `NNN` is counted but not filtered on, so
centromeres that contain assembly gaps are retained: 384 of the 4,551 analysed
chromosomes (8.4%).

### Step 3: the one curated contig

HG03486 hap1 has two `CM` contigs assigned to chr8: `HG03486#1#CM091976.1` (the
chromosome) and `HG03486#1#CM091971.1`, a 90,381 bp contig that is not
chromosome-length. The regions were built on the first; the `.` rows of step 1
name the last, so step 2 writes two extra rows (`subtelo_p` and `centromere`)
for `CM091971.1`, holding only zero counts. Yi-An Chen deleted these two rows by
hand. `aggregate_qc.py --curate` drops them with an explicit filter
(`EXCLUDE`). The chr8 verdict rests on the `CM091976.1` row either way.

## Run

```bash
bash download_prerelease_qc.sh qc_calls          # ~1,386 small BEDs
python3 build_master_qc_table.py --qc-index prerelease_qc_index.tsv --qc-dir qc_calls \
    --contigs inputs/cm_contigs.tsv --regions inputs/regions.tsv -o master_qc_table.full.tsv
python3 aggregate_qc.py master_qc_table.full.tsv > aggregate_qc_v4.tsv
python3 aggregate_qc.py master_qc_table.full.tsv --curate > aggregate_qc_v4.curated.tsv
Rscript ../scripts/QCfilter_apply.R <pangenome.ALLchr.centromere...region.bed> \
    aggregate_qc_v4.curated.tsv centromere Err=0 COLLAPSE=0 COLLAPSE_VAR=0
```

Step 1 takes about 4 s and under 30 MB of memory. The step 4 input BED is the
output of step 1 of `../scripts/HPRC2_centromere_analysis.sh` (6,138 contigs).

## Expected output

Checked on 2026-10-20 from this directory, with the QC BEDs at the listed
sha256:

| file | data rows | lines | sha256 |
|---|---|---|---|
| `master_qc_table.full.tsv` | 225,078 | 225,079 | `2b7faf0a5e43edfac5d4cd296971edb65d92f54c311de99b85f21f600633f218` |
| `aggregate_qc_v4.tsv` | 18,416 | 18,417 | `cccc370f054d644933f69851ac289d117670479b2e92a0765df57175bb73e737` |
| `aggregate_qc_v4.curated.tsv` (`--curate`, LF) | 18,414 | 18,415 | `29ac04bff15bc176f9cd6e2ce293e94484cb08bb539112922f676817bf5ed33f` |
| `aggregate_qc_v4.manual_curation.tsv` (Yi-An Chen; CRLF, no terminator on the last row) | 18,414 | 18,414 | `e149eed58c72cfecc832f9feaf93455fc6b23c66cadcdc636cd41f6cc033f64e` |

- `master_qc_table.full.tsv` is byte-identical to the table the pipeline wrote
  on 2025-11-19.
- `aggregate_qc_v4.tsv` is byte-identical to `gawk -f original/aggregate_qc_v2.awk`
  run on it (with the 4 provenance lines in place or removed).
- The curated table differs from Yi-An Chen's file (the one named in
  `../README.md`) only in line endings: `tr -d '\r'` on that file plus one
  final newline gives `tables/aggregate_qc_v4.curated.tsv` byte for byte. R's
  `read.delim` reads the two identically (`identical()` is `TRUE`).
- Step 4 on either curated file gives the same `*.pass.bed`, 4,650 contigs,
  sha256 `eefde8493f147e01db01723020a56555c5f6d6318a3ac497ce75386f053b6e06`
  (the value in `../README.md`).

## Sensitivity to the QC region

The QC region is not the centromere interval that was clustered (the min/max of
the 100 kb bins of the clustering input), but it is close:

- the two agree to within 1 Mb at both ends for 4,546 of the 4,551 analysed
  chromosomes;
- evaluating QC over the clustering interval instead would remove 4 of the
  4,551 and add 5 chromosomes that currently fail, so the analysed set would go
  from 4,551 to 4,552;
- over all 6,138 chromosome-level centromeres, 11 change status (6 pass to
  fail, 5 fail to pass), so 4,650 would become 4,649.

These counts use the same calls and the same overlap rule as step 1, over the
min/max of each contig's rows in the step 1 BED of
`../scripts/HPRC2_centromere_analysis.sh`. The script that computed them read
the per-region BEDs of the 2025 run and is not included; the same rule over
`inputs/regions.tsv` reproduces the 4,650 exactly.

`sensitivity/` compares the gap-masked Flagger calls used here with the
released ones.

## Provenance

`original/` (each file starts with 4 added comment lines; drop them to get the
original bytes, sha256 given in the header):

- `pangenome.smk`, `get_regions.py`: KaryoScope-archive `workflow/rules/` and
  `workflow/scripts/` at commit `0bf3a7d` ("master pangenome QC table"). Rules
  `get_regions`, `find_errors` (`bedtools intersect -a <QC BED> -b <region BED>`)
  and `aggregate_master_qc_table`. Helper functions such as `get_QC_wildcards`
  are in `common.smk` of that pipeline and are not copied.
- `aggregate_master_qc.py`: the script of `aggregate_master_qc_table`. It was
  never committed; this is the working-tree copy dated 2025-11-19.
- `aggregate_qc_v2.awk`: step 2 as first run (gawk; arrays of arrays). It was
  never committed.

The Snakemake job records of the 2025 run (inputs, commands, input sha256) were
checked against this rebuild: the 99,792 intersection files, the 33,264 region
BEDs, the 1,386 QC BEDs and the 462 `.fai` files all match. Re-running
`find_errors` with bedtools 2.31.1 gives the same 99,792 files.
