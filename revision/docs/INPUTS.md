# Where the inputs live

Every script here names its inputs by absolute path on the TGen cluster. This is
the index. Nothing in this table is in the repository.

## Databases

**The databases are not in `~/.karyoscope/db`.** Always pass
`--db-root $KS_SOFTWARE/github/barthel/KaryoScope/resources/registry`.

| database | what it is | used by |
|---|---|---|
| `KS_human_CHM13` | **v1**, KMC backend | the published centromere-clustering and D4Z4 inputs, and Supp. Fig. 5 panel E |
| `KS_human_CHM13_v2` | v2, KMC backend | the concordance work, A4, A5 |
| `HKS_human_CHM13_v2` | v2, HKS backend | the runtime benchmarks |
| `HKS_human_CHM13_v2_k{11,21,31,41,51,61}` | the *k* sweep rebuilds (**v2** vocabulary) | `a2-k-sensitivity` (under `$KS_WORK/ks_build_sources/chm13v2_ksweep/db`) |
| `HKS_human_CHM13_v1like_region_k{21,31,51}` | region-only rebuilds of the **v1** feature set at other *k* | `a2-k-sensitivity/step4` — the D4Z4 and centromere analyses (under `$KS_WORK/ks_build_sources/v1like_ksweep/db`). Built with `build: external_memory:` (20 GB peak, not ~75 GB). The *k* = 31 member reproduces the installed v1 KMC database bit-for-bit on the D4Z4 windows and reproduces the published centromere clustering on 4,551/4,551 haplotypes. |

**v1 and v2 are not interchangeable, and this bites anything that reuses the six-*k*
family.** The `region` vocabularies genuinely differ (v1 uses `*_specific` suffixes,
`hor_multigroup1`), so a v2 annotation cannot be relabelled into v1's namespace.
Worse, the *intervals* differ too: v2 labels the whole 117 kb D4Z4 array as one
`bSat` interval where v1 leaves it `q_arm`, so the published D4Z4 unit count
collapses 28,580 -> 1,378 under any v2 database, and the centromere clustering at
*k* = 31 gives 1,195 minor haplotypes against the published 821. Anything
reproducing a published centromere or D4Z4 result needs v1 — regenerate the v1
region set at other *k* per `a2-k-sensitivity/step4/RESULTS.md`.

## Assemblies

| what | where |
|---|---|
| HG002 v1.1, both haplotypes and combined diploid | `$KS_PROJECTS/` — see the per-script paths |
| T2T-CHM13v2.0 | as above |
| the three Robertsonian samples (GM03417, GM03786, GM04890) | `$KS_GUARRACINO/ks-robs/results/<sample>/assembly/haplotype{1,2}/` |
| HPRC Release 2 pangenome | per-chromosome extracts, below |

## Pangenome annotation extracts

One directory carries the inputs for both the centromere-clustering and D4Z4
analyses, annotated against **v1**:

```
$KS_PROJECTS/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed/
```

* `pangenome.chr{1..22,X,Y}.centromere.KS_human_CHM13.presmoothed.region.bed` —
  24 files, the centromere-clustering input. Contig column is the **bare**
  `HG00097#1#CM094060.1` form.
* `pangenome.chr{4,10}.subtelo_q.KS_human_CHM13.presmoothed.{region,subtelomeric,gene,repeat}.bed` —
  8 files, the D4Z4 input. Contig column is the **prefixed**
  `chr4_hap2_HG00099#2#CM087359.1` form.

The two naming conventions are a real difference between the pipelines, not a
fault in either. `QCfilter_apply.R` joins on contig name as a plain string, so a
mismatch yields an empty output rather than an error — check row counts.

## Other inputs

| what | where | note |
|---|---|---|
| NucFlag QC curation | `$KS_WORK/aggregate_qc_v4.manual_curation (2).tsv` | 18,413 rows. **Manual curation — cannot be regenerated.** |
| satellite colour palette | `resources/registry/KS_human_CHM13/colors.tsv` | the driver's default `resources/databases/...` path is the old layout |
| TAR1 / ITS element BEDs | S3, sha256 `abd3197…` | **not** the UCSC bigBed the `repeat` feature set uses |
| HMM D4Z4 counts for chr4 (Fig. S8B) | `barthel-lab/HPRC-D4Z4:Hailey-notes/HPRC_chr4_counts.xlsx` | 534 contigs; joined on the contig id it reproduces the published *r* = 0.998 at *n* = 231. Not vendored here. The derived per-*k* correlation table is `a2-k-sensitivity/step4/results/d4z4_hmm_correlation_by_k.tsv`, and `scripts/job3-d4z4/run_hmm_counts.sh` re-derives HMM counts independently from the Methods recipe (agrees at *r* = 0.9991) |
| the 582 D4Z4 and 6,138 centromere window coordinates | `a2-k-sensitivity/step4/results/{d4z4,cen}_windows.tsv` | extracted from the release 2 FASTAs with `samtools faidx`; every published window is forward-frame because the archive's window scripts silently drop reverse-complemented contigs |
| the *k* = 31 scaffolded chromosome BEDs | `$KS_WORK/a2_step3/sweep/<sample>/k31/` | the ROB sweep reads these, so it needs no annotation step |

## Sample naming

The manuscript names two validated cell lines by **Coriell catalogue number**
(`GM21144`, `GM18960`); the assemblies use the **1000G/HPRC name** (`NA21144`,
`NA18960`). The pangenome tables hold `HG` and `NA` samples and **no `GM` at
all**, so a lookup by the manuscript's name returns "absent" rather than an
error.

## Toolchain

| tool | where | note |
|---|---|---|
| `karyoscope` and the analysis env | `$KS_SOFTWARE/miniforge3_new/envs/karyoscope/bin` | use the absolute path in SLURM scripts |
| `hks` | `$KS_SOFTWARE/github/barthel/HKS/target/release/hks` | via `KARYOSCOPE_HKS` |
| `Rscript` | `$KS_HOME/software/bin/Rscript` | needed by the centromere pipeline's steps 1–2 |
| LaTeX | `$KS_WORK/containers/texlive.sif` | singularity; there is no `pdflatex` on the host |

## A7 review checks (2026-09-07)

Nothing new was annotated except three contigs and their reverse complements, and one gene-only database.

| input | path | job |
|---|---|---|
| Ey15-2 and Col-CEN nuclear A3 annotations, **both stages kept** | `$KS_WORK/a3_crossgenome/annot_{ey15_2,colcen_nuclear}/` (job 35378568) | 1 |
| HG002 hap1 / CHM13 A4 annotations, all six sets, both stages | `$KS_WORK/a4_smoothing/annot_{hg002hap1,chm13}/` (job 35378688); N-mask `a4_smoothing/HG002.hap1.n_regions.renamed.bed` | 2, 4 |
| recipe's gene build input | `$KS_WORK/recipe_repro/chm13_gene.bed` (sha256 03ae30a1…, 612,909 rows) | 2 |
| A4 `hks node-stats` output | `$KS_WORK/a4_smoothing/raw_nodestats/nodestats_shipped_gene.txt` | 2 |
| gene-only database **without** priorities | `$KS_WORK/a7_review_checks/job2_noprio/db/HKS_human_CHM13_v2_gene_noprio` | 2 |
| Job 3 query contigs, both orientations, and their HKS/KMC annotations | `$KS_WORK/a7_review_checks/job3/{fa,annot_hks,annot_kmc}/` | 3 |
| frozen k = 31 D4Z4 caller table (reproduces the census bit-for-bit) | `$KS_WORK/a2_step4/job3_d4z4/base_published/results/d4z4_array_analysis.tsv`; windows `a2_step4/common/d4z4_windows.tsv` | 5 |
