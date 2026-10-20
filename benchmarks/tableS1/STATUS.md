# Table S1 cold/warm re-run — status / pick-up notes

**Job 38434006** (`ks_tableS1_cold`, clong, 48 h limit), submitted 2026-10-13, **STARTED 14:28 on g-c-1-4-34**, ETA ~00:30-01:30 Oct 14. First run (diploid HKS r2 cold) checked: 6 feature sets, exit 0, resident 0, others=0, extractor parses it.
Expected runtime ~9-10 h (old Table S1 job: 18 runs in 4h11; now 36 runs + eviction).
Script `run_tableS1.sh` (design in its header): same protocol as `../threadsweep_cold`
(one node, 42 CPU/490 GB, random order seed 20261013, evict + fincore-verify, cold then warm),
Table S1 design (3 inputs x 2 backends x 3 reps, all six feature sets, -t 16), pinned code from
`../threadsweep_cold/` (ks-v2.3.0, hks-0.4.0, get_featureIDs). Resumable.

Requested by the user after the other (manuscript) Claude asked for it; its asks:
cold+warm for every Table S1 column, exact input sizes, the node name.

## Facts already settled (2026-10-13)
- Old Table S1 job 35247555 ran on **g-c-1-7-24**. g-c-1-4-02 was the CANCELLED 64G attempt 35247205.
- Old run order was rep-major; rep 1 was slowest for all 6 input x backend (8-20 %): mixed cache state.
  Old per-run walls survive only in a transcript (results_35247555.tsv), the logs were purged.
- Input sizes (.fai): CHM13 v2.0 3,117,292,070 bp (25 seqs incl. chrM);
  HG002 v1.1 hap1 2,947,895,164 bp (23 seqs) -> **2.95 Gb, Table S1 caption's 2.89 is wrong**;
  HG002 v1.1 diploid 5,999,424,718 bp (47 seqs incl. chrM).
- Input files: CHM13 = uncompressed `chm13v2.0.fasta` (the old script's choice is unrecoverable);
  HG002 = the bgzipped .fasta.gz files.
- Outputs go to node-local /tmp (old run's output location unknown; compression column may differ).

## When it finishes
1. `tail tableS1_38434006.log` -> ALL DONE; `wc -l runs/index.tsv` = 37; residency 0 on cold rows; others=0.
2. `$KS_WORK/d4z4_qc/env/bin/python extract_tableS1.py` -> tableS1_runs.tsv,
   tableS1_summary.tsv (cold + warm; check feature_sets == 6 and exit == 0 in tableS1_runs.tsv).
3. Add the results to `../threadsweep_cold/HANDOFF-figS6.md` (or a Table S1 hand-off) for the manuscript Claude.

## Progress check 2026-10-13 16:05
6/36 runs done (points 0-2), all rc=0, resident 0 before every cold run, others=0. ETA unchanged (~00:30-01:30 Oct 14).
Interim (n=1, min): diploid HKS cold 15.17 / warm 13.23; diploid KMC cold 23.74 / warm 21.40;
CHM13 HKS cold 7.31 / warm 5.25. Cold penalty sits in assignment only (~2 min); smoothing,
compression and RSS identical cold vs warm, as in the Fig. S6 sweep.

## TOMORROW: pick-up checklist
1. `sacct -j 38434006 -X -o State,Elapsed`; `tail -2 tableS1_38434006.log` -> ALL DONE.
   If it died: `sbatch run_tableS1.sh` resumes (skips points with a warm row).
2. `wc -l runs/index.tsv` = 37; `grep -o 'others=[0-9]*' tableS1_*.log | sort | uniq -c` all 0;
   cold rows resident_bytes_before = 0.
3. `$KS_WORK/d4z4_qc/env/bin/python extract_tableS1.py`; check tableS1_runs.tsv
   feature_sets == 6 and exit == 0 everywhere.
4. Write the results (warm = Table S1 replacement; cold = first-run cost for a full annotation;
   old-vs-new comparison; old table values are in the hand-off/memory) into the hand-off for the
   manuscript Claude: `../threadsweep_cold/HANDOFF-figS6.md` §5 item 6 + §6 rows, or a new
   `HANDOFF-tableS1.md` next to this file. Also recompute the main-text numbers that cite Table S1
   (main.tex ~196: 13.7 +- 0.4 / 24.8 +- 1.8 KMC, 7.9 +- 0.8 / 15.2 +- 1.1 HKS, 10.1 / 34.9 / 43.7 GiB;
   Supp Note 2 ~676: 1.73x hap / 1.63x diploid / 3.4x memory) and Table S2's KaryoScope rows
   (7.9 / 13.7 min, "n=3 means from Table S1").
5. Nothing is committed or pushed; the paper checkout is SHARED (use a git worktree if editing).

## Progress check 2026-10-13 21:39
34/36 done, all clean. Cold CHM13 HKS r2 outlier (assign 949 s vs ~388 s) = shared-FS I/O stall: same user CPU + same bytes read, 403% vs ~930% CPU. Warm unaffected.

## DONE 2026-10-14
Job finished 22:00:41 Oct 13: ALL DONE, 36/36 rc=0, 6 feature sets each, others=0, cold resident 0.
Extracted; hand-off written: `HANDOFF-tableS1.md` (warm table + LaTeX rows, cold disclosure, old-vs-new,
claim-by-claim updates vs origin/revision/fpb-comments-2026-10-10). Fig S6 hand-off now points to it.
