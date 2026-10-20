# Table S1 re-run: results and hand-off (2026-10-14)

For the Claude instance working on the KaryoScope manuscript. Supplementary Table S1
(`\stref{tab:S1_runtime_memory}`, runtime and memory of `karyoscope annotate`) was re-run
under the same cache protocol as the new Fig. S6 sweep
(`$KS_WORK/ks_slurm/threadsweep_cold/HANDOFF-figS6.md`). The user chose
**warm only** for Fig. S6, and Table S1 is reported the same way: **the warm numbers
replace the table.** The cold numbers are a one-sentence disclosure.

Nothing has been committed or pushed. The paper checkout at
`$KS_SOFTWARE/github/barthel/KaryoScope-paper` is SHARED, so make edits in a
`git worktree`. Line numbers below are from `origin/revision/fpb-comments-2026-10-10`
(`manuscript/main.tex`), which already contains the Fig. S6 rewrite.

---

## 1. Why it was re-run

The original Table S1 (job 35247555, node g-c-1-7-24) ran its 18 runs back-to-back with no
cache control. Replicate 1 was the slowest run for every input × backend (8–20%), so the
table mixed first-run and repeat-run timings. The original per-run logs were purged.

## 2. What was run

| Item | Value |
|---|---|
| Job | 38434006, one job on one node (**g-c-1-4-34**), 2026-10-13 14:28 to 22:00. All 36 runs exited 0 |
| Hardware | 2× Intel Xeon Gold 6342 (Ice Lake, 2.80 GHz), 48 cores, SMT off, 503 GiB. Same node class as the original Table S1 and the new Fig. S6 sweep. 42 CPUs and 490 GB were allocated; **no other job shared the node during any run** (`others=0` for all 36) |
| Software | KaryoScope v2.3.0 and hks 0.4.0, the same pinned copies as the Fig. S6 sweep |
| Databases | `KS_human_CHM13_v2` (KMC), `HKS_human_CHM13_v2` (HKS) |
| Inputs | T2T-CHM13v2.0 (3,117,292,070 bp); HG002 v1.1 hap1 (2,947,895,164 bp); HG002 v1.1 diploid (5,999,424,718 bp) |
| Design | 3 inputs × 2 backends × 3 replicates = 18 points, all six feature sets, 16 threads, run in a fixed random order (seed 20261013) |
| Per point | Evict the database, input, binaries and source from the OS page cache and confirm with `fincore` that nothing is resident (it was 0 bytes before every cold run). Then a **cold** run, then an identical **warm** run immediately afterwards |
| Outputs | Written to node-local `/tmp` and deleted between runs. Every run wrote all six feature sets |

## 3. Results: WARM (the replacement Table S1)

Minutes, mean ± s.d., n = 3. Peak RSS in GiB.

| Assembly | Backend | Assignment | Smoothing | Compression | End-to-end | Peak RSS |
|---|---|---|---|---|---|---|
| T2T-CHM13v2.0 | HKS | 4.34 ± 0.04 | 0.60 ± 0.01 | 0.29 ± 0.01 | **5.24 ± 0.03** | 10.15 ± 0.01 |
| T2T-CHM13v2.0 | KMC | 3.28 ± 0.04 | 8.06 ± 0.05 | 0.29 ± 0.00 | **11.64 ± 0.02** | 34.42 ± 0.15 |
| HG002 hap1 | HKS | 5.28 ± 0.03 | 0.67 ± 0.00 | 0.30 ± 0.00 | **6.25 ± 0.03** | 10.07 ± 0.00 |
| HG002 hap1 | KMC | 2.88 ± 0.23 | 8.27 ± 0.06 | 0.32 ± 0.02 | **11.48 ± 0.16** | 34.85 ± 0.12 |
| HG002 diploid | HKS | 10.67 ± 0.36 | 1.46 ± 0.00 | 0.68 ± 0.01 | **12.82 ± 0.36** | 10.10 ± 0.01 |
| HG002 diploid | KMC | 5.33 ± 0.16 | 15.05 ± 0.04 | 0.86 ± 0.04 | **21.27 ± 0.13** | 43.81 ± 0.14 |

The LaTeX rows, ready to paste over the table body (same column order as the current table):

```latex
T2T-CHM13v2.0 & HKS & $4.34 \pm 0.04$ & $0.60 \pm 0.01$ & $0.29 \pm 0.01$ & $5.24 \pm 0.03$ & $10.15 \pm 0.01$ \\
T2T-CHM13v2.0 & KMC & $3.28 \pm 0.04$ & $8.06 \pm 0.05$ & $0.29 \pm 0.00$ & $11.64 \pm 0.02$ & $34.42 \pm 0.15$ \\
HG002 hap1 & HKS & $5.28 \pm 0.03$ & $0.67 \pm 0.00$ & $0.30 \pm 0.00$ & $6.25 \pm 0.03$ & $10.07 \pm 0.00$ \\
HG002 hap1 & KMC & $2.88 \pm 0.23$ & $8.27 \pm 0.06$ & $0.32 \pm 0.02$ & $11.48 \pm 0.16$ & $34.85 \pm 0.12$ \\
HG002 diploid & HKS & $10.67 \pm 0.36$ & $1.46 \pm 0.00$ & $0.68 \pm 0.01$ & $12.82 \pm 0.36$ & $10.10 \pm 0.01$ \\
HG002 diploid & KMC & $5.33 \pm 0.16$ & $15.05 \pm 0.04$ & $0.86 \pm 0.04$ & $21.27 \pm 0.13$ & $43.81 \pm 0.14$ \\
```

The "0.00" SDs are real: below 0.005 min (0.3 s).

**Derived ratios (warm, end-to-end and peak RSS):**

| Input | HKS speed-up over KMC | KMC/HKS memory |
|---|---|---|
| T2T-CHM13v2.0 | 2.22× | 3.39× |
| HG002 hap1 | 1.84× | 3.46× |
| HG002 diploid | 1.66× | 4.34× |

So: **HKS is 1.7–2.2× faster at 3.4–4.3× less peak memory** (the current text says
1.6–1.7× and 3.4–4.3×).

## 4. Results: COLD (not for the table; one sentence of disclosure)

| Assembly | Backend | Assignment | End-to-end | Peak RSS |
|---|---|---|---|---|
| T2T-CHM13v2.0 | HKS | 9.59 ± 5.39 (median 6.54) | 10.56 ± 5.52 (median 7.43) | 10.12 |
| T2T-CHM13v2.0 | KMC | 4.75 ± 0.66 | 13.18 ± 0.59 | 34.44 |
| HG002 hap1 | HKS | 6.84 ± 0.48 | 7.82 ± 0.48 | 10.08 |
| HG002 hap1 | KMC | 3.60 ± 0.29 | 12.20 ± 0.29 | 34.86 |
| HG002 diploid | HKS | 12.19 ± 0.72 | 14.37 ± 0.71 | 10.11 |
| HG002 diploid | KMC | 7.42 ± 0.33 | 23.55 ± 0.43 | 43.65 |

- **First-run cost:** the median extra assignment time per point (cold minus warm) is
  **88–130 s for HKS and 38–118 s for KMC**. That is larger than in the Fig. S6 sweep
  (41 / 57 s) because Table S1 loads all six feature sets (the HKS index is 25.8 GB on
  disk) rather than one.
- **Smoothing, compression and peak RSS are the same cold and warm.** Only assignment changes.
- **One cold outlier, CHM13 HKS replicate 2:** assignment took 949 s, against 385–392 s for
  the other two replicates. It used the same user CPU time (about 4,057 s) and read the same
  25.8 GB, but averaged 403% CPU against about 930%, so it was waiting on the shared network
  filesystem. It is not a code effect, and the warm run that followed was normal (258 s).
  If cold numbers are ever tabulated, report medians or footnote this run.

Suggested disclosure, matching the Fig. S6 paragraph (line 670): *Each configuration was
run twice in succession after clearing the database, input and software from the
operating-system file cache, and the second run is reported. On the first run, loading the
six-feature-set database added a median of 88–130 s (HKS) and 38–118 s (KMC) to assignment;
smoothing, compression and peak memory were unaffected.*

## 5. Old vs new (warm), for the record

| | Old end-to-end | New (warm) | Change |
|---|---|---|---|
| CHM13 HKS | 8.67 ± 0.33 | 5.24 ± 0.03 | −40% |
| CHM13 KMC | 13.57 ± 0.64 | 11.64 ± 0.02 | −14% |
| hap1 HKS | 7.91 ± 0.82 | 6.25 ± 0.03 | −21% |
| hap1 KMC | 13.71 ± 0.44 | 11.48 ± 0.16 | −16% |
| diploid HKS | 15.23 ± 1.05 | 12.82 ± 0.36 | −16% |
| diploid KMC | 24.77 ± 1.80 | 21.27 ± 0.13 | −14% |

Peak RSS is unchanged within 0.1 GiB everywhere. The drop has two sources:
1. **Assignment:** the old runs mixed cold and warm cache states.
2. **Compression, and part of smoothing:** these were also faster and far less variable
   (compression 0.49–2.12 min with SDs up to 1.19 old, against 0.29–0.86 min with SDs ≤ 0.04
   now). The new runs wrote to node-local disk. Where the old runs wrote is unrecorded, but
   their spread is consistent with shared network storage. Both phases write output, and
   the cache protocol cannot change them.

So the new numbers measure computation on the stated node with local output. Saying that
outputs went to node-local disk keeps the table honest.

## 6. Claim-by-claim updates

| Location | Current | New (warm) |
|---|---|---|
| **Main text 198** | KMC hap1 $13.7 \pm 0.4$ min, $34.9 \pm 0.2$ GiB | **$11.5 \pm 0.2$ min**, $34.9 \pm 0.1$ GiB |
| Main text 198 | KMC diploid $24.8 \pm 1.8$ min at $43.7 \pm 0.1$ GiB | **$21.3 \pm 0.1$ min at $43.8 \pm 0.1$ GiB** |
| Main text 198 | HKS $7.9 \pm 0.8$ min per haplotype, $15.2 \pm 1.1$ diploid, 10.1 GiB | **$6.3 \pm 0.03$** (or $6.25 \pm 0.03$) **min**, **$12.8 \pm 0.4$ min**, 10.1 GiB |
| Table S1 body (738–745) | old rows | §3 LaTeX rows |
| **Table S1 caption 734** | "HG002~v1.1 haplotype~1, 2.89~Gbp" | **2.95~Gbp** (`.fai`: 2,947,895,164 bp; Table S2 caption and Methods 472 already say 2.95 / 2,947,895,164) |
| Table S1 caption 734 | no cache condition stated | add: each configuration run twice in succession after clearing the OS file cache, second run reported; outputs written to node-local disk |
| **Table S2 rows 760–761** | HKS 7.9 min / 10.1; KMC 13.7 min / 34.9 | **HKS 6.3 min / 10.1; KMC 11.5 min / 34.8** (34.85 rounds to 34.8 or 34.9; main text uses 34.9, so keep the two consistent) |
| Supp Note 2, 666 | "executed sequentially within a single job" | still true; add the cache protocol and the random order (or point to §4 wording above) |
| **Supp Note 2, 668** | hap1 $13.71 \pm 0.44$ KMC / $7.91 \pm 0.82$ HKS at $34.85 \pm 0.16$ / $10.12 \pm 0.01$ | **$11.48 \pm 0.16$ / $6.25 \pm 0.03$** at **$34.85 \pm 0.12$ / $10.07 \pm 0.00$** |
| Supp Note 2, 668 | diploid $24.77 \pm 1.80$ / $15.23 \pm 1.05$ at $43.69 \pm 0.11$ / $10.19 \pm 0.05$ | **$21.27 \pm 0.13$ / $12.82 \pm 0.36$** at **$43.81 \pm 0.14$ / $10.10 \pm 0.01$** |
| Supp Note 2, 668 | CHM13 $13.57 \pm 0.64$ / $8.67 \pm 0.33$ at $34.45 \pm 0.08$ / $10.10 \pm 0.04$ | **$11.64 \pm 0.02$ / $5.24 \pm 0.03$** at **$34.42 \pm 0.15$ / $10.15 \pm 0.01$** |
| Supp Note 2, 668 | "1.6--1.7$\times$ faster … 3.4--4.3$\times$ less peak memory" | **1.7--2.2$\times$** faster; 3.4--4.3$\times$ unchanged |
| Supp Note 2, 668 | RepeatMasker "approximately 190$\times$ longer than KaryoScope-KMC" | **approximately 230$\times$** (2,634 min / 11.48 min = 229×; against HKS it is 421×) |
| Supp Note 2, 672 (last sentence) | HKS assignment faster for repeat alone, slower for six sets | still true (hap1 assignment: HKS 5.28 vs KMC 2.88 min) |

If origin/main is used instead, its old version of Supp Note 2 also has "under fourteen
minutes … under eight minutes" (main ~674) and "1.73× / 1.63× / 3.4×" (main ~676). Under the
new numbers these become **"under twelve minutes … under seven minutes"** and
**1.84× / 1.66× / 3.5×** (hap1 memory ratio 3.46). The revision branch has already rewritten
that passage.

## 7. File locations

All under `$KS_WORK/ks_slurm/tableS1_cold/`:

| File | What it is |
|---|---|
| `run_tableS1.sh` | sbatch script, with the design rationale in its header |
| `manifest.tsv` | the 18 points in their random run order |
| `runs/index.tsv` | one row per run (36): input, backend, rep, mode, node, job, bytes resident before the run, start time, load average, log path |
| `runs/*.log` | full `karyoscope -v` and `/usr/bin/time -v` log for each run |
| `tableS1_38434006.log` | job log, including the per-run co-tenant count (`others=`) and `lscpu` |
| `extract_tableS1.py` | log-to-TSV extractor |
| **`tableS1_runs.tsv`** | per-run values: every phase, end-to-end, RSS, feature-set count (6 in every row), exit code (0 in every row) |
| **`tableS1_summary.tsv`** | means and SDs (minutes and GiB), cold and warm |
| `STATUS.md` | working notes and history |
