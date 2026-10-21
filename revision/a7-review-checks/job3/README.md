# A7 Job 3: reverse-complement symmetry of smoothing

`job3_hks.sh` extracts the query sequences, reverse-complements them with
`seqtk seq -r`, annotates both orientations in one `karyoscope annotate` run
(`HKS_human_CHM13_v2`, all six feature sets, `hks smooth`), and runs
`job3_compare.py`, which maps each reverse-complement interval back to forward
*k*-mer coordinates and counts differing bases per feature set and stage.
`job3_kmc.sh` repeats the annotation with the KMC database `KS_human_CHM13_v2`
(smoothing in `karyoscope/core/smooth.py`). Results go to
`$KS_WORK/a7_review_checks/results/rc_symmetry_{hks,kmc}*.tsv`; the RESULTS
write-up is `../RESULTS.md`, Job 3.

## What it backs in the manuscript

The Methods sentence on smoothing: "We tested this by annotating three
sequences (a 44 Mb acrocentric chromosome and a 67 Mb autosome of HG002
haplotype 1, and the 133 Mb Robertsonian fusion scaffold of GM03417) and their
reverse complements. The labellings were identical in all six feature sets on
both backends."

## Public and controlled parts

- **HG002 (public).** With `ROBS_CONFIG` unset, the scripts run the two HG002
  haplotype 1 sequences (`chr21_hap1`, `chr20_hap1`) only. The result tables in
  `../results/rc_symmetry_*.tsv` are this HG002 part.
- **GM03417 (dbGaP phs003920).** With `ROBS_CONFIG` set to a filled ROB
  configuration, the `a7_job3` section adds the fusion contig (sample,
  haplotype, contig name in the original assembly, file stem, description) as
  the third sequence, read from `$KS_GUARRACINO/ks-robs/results/`.

## Run

```bash
export KS_WORK=... KS_SOFTWARE=... KS_REFERENCES=...   # + KS_GUARRACINO and ROBS_CONFIG for GM03417
mkdir -p $KS_WORK/a7_review_checks/{scripts,results,logs}
cp job3_compare.py $KS_WORK/a7_review_checks/scripts/
jid=$(sbatch --parsable job3_hks.sh)          # about 3 min, 10 GB
sbatch --dependency=afterok:$jid job3_kmc.sh  # about 25 min, 19 GB
```

## Changes from KaryoScope-revisions `b7d2092`

The fusion contig is read from `$ROBS_CONFIG` instead of being written in the
scripts, and the ROB sequence is optional; cluster paths are replaced by the
variables in the repository README (Paths). Checked 2026-10-20 (see the
repository README, Controlled-access data).

## Checked 2026-10-20

The scripts here read the GM03417 contig from the ROB configuration instead of
naming it. Run with a filled configuration, they write all six
`rc_symmetry_*.tsv` tables byte-identical to the original scripts' output, which
is in turn byte-identical to the tables behind the Methods sentence. With
`ROBS_CONFIG` unset, the six tables are byte-identical to the HG002-only tables
in `../results/`.
