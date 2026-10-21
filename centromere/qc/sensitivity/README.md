# Gap-masked vs released Flagger calls (centromere regions)

The QC table in `../tables/` was built from the HMM-Flagger v1.2.0 conservative
calls in their per-haplotype, gap-masked form
(`*.conservative.canonical.no_Hap.hap{1,2}.bed`). The released product on the
same S3 path is `*.conservative.bed`. These two scripts compare them over the
same centromere regions (`../inputs/regions.tsv`), with the same aggregation
and the same filter (`Err=0 COLLAPSE=0 COLLAPSE_VAR=0`).

| script | output |
|---|---|
| `compare_prerelease_released.py` | `comparison.tsv` (calls and bp per tool and flag, both call sets), `passfail_diff.tsv` (chromosomes whose verdict differs), summary on stderr |
| `addendum_gaps_outliers.py` | `addendum.txt`, `addendum_regions.tsv`, `addendum_nnn_per_chrom.tsv`, `addendum_gap_vs_outlier.tsv`: gap-containing centromeres among the 4,551 analysed, gap vs Outlier (Fisher tests), and how much of each gap the released `Err` calls cover |

## Inputs

- the pre-release QC BEDs (`../download_prerelease_qc.sh <dir>`, then `--qc-dir <dir>`);
- the released Flagger BEDs: `download_released_flagger.sh <dir>` fetches the
  462 files of `released_flagger_index.tsv` (HiFi and ONT, one file per sample
  holding both haplotypes) and checks their sha256. The released NucFlag BEDs
  are byte-identical to the pre-release ones, so both call sets use the same
  NucFlag files;
- `../tables/aggregate_qc_v4.curated.tsv` (`--curated`), used to check that the
  pre-release recount equals every centromere row of the table;
- from `../../scripts/HPRC2_centromere_analysis.sh`: the step 1 BED
  (`--clustering-bed`), step 3 `allchr_structure.sequence_assignments.tsv`
  (`--assignments`), step 5 `allchr_outlier_annotations.tsv` (`--gated`) and
  step 7b `allchr_allele_pairs.tsv` (`--allele-pairs`);
- `--clustered`: the step 3 assignments of the 4,551 analysed chromosomes, as a
  TSV with columns `sample haplotype chrom contig cluster cluster_type` (step 3
  output joined to the QC-passed contigs, minus chr13, chr14, chr15, chr21,
  chr22 and chrY);
- `--gaps`: a BED of the assembly N runs on those contigs, for example
  `seqtk gap -l 1 <assembly.fa.gz>` restricted to the contigs.

```bash
python3 compare_prerelease_released.py --qc-index ../prerelease_qc_index.tsv --qc-dir <pre> \
    --released-dir <rel> --regions ../inputs/regions.tsv --curated ../tables/aggregate_qc_v4.curated.tsv \
    --clustering-bed <step 1 BED> --out-dir <out>
python3 addendum_gaps_outliers.py --qc-index ../prerelease_qc_index.tsv --qc-dir <pre> --released-dir <rel> \
    --regions ../inputs/regions.tsv --clustered <clustered.tsv> --gaps <gaps.bed> \
    --assignments <...sequence_assignments.tsv> --allele-pairs <allchr_allele_pairs.tsv> \
    --gated <allchr_outlier_annotations.tsv> --out-dir <out>
```

The second script needs pandas and scipy (`../../../environment/centromere.yml`).

## Result (2026-10-20)

- The pre-release recount equals all 6,138 centromere rows of the QC table and
  keeps 4,650.
- The released calls keep 4,596: 54 chromosomes change, all from pass to fail,
  52 of them among the 4,551 analysed.
- In all 54 the only difference is assembly gaps: the gap-masked files label
  `Err` that falls in a gap as `NNN`, which is not filtered on; the released
  files keep it as `Err`. NucFlag `COLLAPSE` and `COLLAPSE_VAR` never change.
- 384 of the 4,551 analysed centromeres contain an assembly gap (NNN > 0), the
  same 384 that contain an N run. They are not enriched among Outlier or minor
  clusters (ungated stage 3: odds ratio 1.01, p = 0.92).
- The released `Err` calls cover gaps of about 16 kb and longer. Shorter gaps
  and 100 bp placeholder gaps fall inside `Hap` calls and carry no flag.

These scripts were re-run from this directory on 2026-10-20; all six output
files and the stderr summary are byte-identical to the first run.
