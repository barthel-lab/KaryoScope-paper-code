# Fig. S6 thread sweep: results and hand-off (2026-10-13)

For a Claude instance working on the KaryoScope manuscript. The sweep behind
Supplementary Fig. S6 ("Thread scaling of the two principal KaryoScope annotation
pipeline steps, by backend"; `\sfref{fig:scaling-sup}`) was re-run from scratch.
**Decision already made by the user: the manuscript shows the WARM data only.**
The text describing the figure (caption, Supplementary Note 2, plot-script
docstring) is now out of date and must be updated to match.

Nothing has been committed or pushed. The paper checkout at
`$KS_SOFTWARE/github/barthel/KaryoScope-paper` is SHARED, so make
edits in a `git worktree`.

---

## 1. Why it was re-run

The published figure came from jobs 35376001/35376043 (`$KS_WORK/ks_slurm/threadsweep/`).
Cache state was confounded with thread count in that sweep:

- **Order:** thread counts ran in descending order, so the one run that read the index from disk landed on t=32. That is the 104 s replicate behind the wide 32-thread HKS error bar the caption explains as "one slow replicate".
- **Second job:** replicate 3 came from a second job and was 5–30% faster at every point.
- **Shared node:** both jobs shared their node with other tenants.

The user asked for every point to be measured the same way, with cache ruled out as
a confounder. The old Supplementary Note paragraph about replicate variation
(main.tex ~line 678) describes that confounded sweep and no longer applies.

## 2. What was run (new sweep)

| Item | Value |
|---|---|
| Job | 38373939, one job on one node (g-c-1-6-21), 2026-10-12 15:54 to 2026-10-13 04:12 |
| Hardware | Intel Xeon Gold 6342 (Ice Lake, 2.80 GHz), the same node class as Table S1. 42 CPUs and 490 GB allocated, the most a single job can get here. **No other job shared the node during any of the 72 runs** (`others=0` logged per run). |
| Software | KaryoScope v2.3.0 (pinned git worktree) and hks 0.4.0 (pinned binary copy) |
| Input | HG002 v1.1 hap1 (2.95 Gb), `repeat` feature set |
| Databases | `HKS_human_CHM13_v2` (HKS backend) and `KS_human_CHM13_v2` (KMC backend). **The current caption names only KS_… for both; that is wrong.** |
| Design | 2 backends × 6 thread counts (1, 2, 4, 8, 16, 32) × 3 replicates = 36 points, run in a fixed random order (seed 20261012) so drift cannot line up with thread count. The 24-thread point is gone (it was already dropped from the plotted figure). |
| Per point | (a) evict the database, input, binaries and source from the OS page cache (`posix_fadvise DONTNEED`) and confirm with `fincore` that nothing is resident; the job aborts otherwise. (b) **cold** run. (c) an identical **warm** run immediately afterwards, with the cache as the cold run left it. |
| Outputs | Written to node-local /tmp and deleted between runs. All 72 runs exited 0. |

"Warm" therefore means: **the database files were already in the operating-system
file cache because the same command had just run.** This is the state of every run
after the first one on a node, e.g. when annotating many assemblies in sequence.
"Cold" means the database (24 GB HKS / 19 GB KMC on disk) is read from network
storage first, as on a user's very first run.

## 3. Results: WARM (the version going in the manuscript)

Mean ± SD, n = 3. Throughput = 2,950 Mb / assignment wall-clock.

| Threads | HKS assign (Mb/s) | KMC assign (Mb/s) | HKS/KMC | Rust smooth (s) | Python smooth (s) | HKS RSS (GiB) | KMC RSS (GiB) |
|---|---|---|---|---|---|---|---|
| 1 | 6.49 ± 0.07 | 2.33 ± 0.01 | 2.79× | 42.2 ± 0.0 | 1302.3 ± 17.9 | 9.36 | 21.66 |
| 2 | 12.45 ± 0.34 | 4.22 ± 0.29 | 2.95× | 25.1 ± 0.2 | 664.0 ± 5.6 | 9.41 | 23.26 |
| 4 | 19.64 ± 0.80 | 7.29 ± 0.16 | 2.69× | 16.6 ± 0.2 | 363.7 ± 4.0 | 9.53 | 26.01 |
| 8 | 33.51 ± 0.16 | 12.41 ± 0.21 | 2.70× | 13.8 ± 0.1 | 200.7 ± 0.6 | 9.71 | 30.28 |
| 16 | 54.83 ± 0.27 | 16.01 ± 0.24 | **3.42×** | 13.8 ± 0.1 | 142.7 ± 1.5 | 10.07 | 34.90 |
| 32 | 79.80 ± 0.54 | 15.50 ± 0.80 | 5.15× | 13.8 ± 0.0 | 126.0 ± 2.6 | 10.83 | 36.08 |

End-to-end wall-clock of one `annotate` run (`repeat` set only, 16 threads, warm):
HKS 71.2 s (1.19 min), KMC 331.5 s (5.53 min), **4.65×**.

**What the warm data show:**

- **HKS assignment** keeps scaling to 32 threads (54.8 → 79.8 Mb/s). **KMC assignment** stops improving after 16 threads (16.0 → 15.5).
- **Error bars are now tiny everywhere.** The largest assignment SD is 0.8 Mb/s, so the "wide interval at 32 threads" no longer exists.
- **Rust smoothing** falls from 42.2 to 13.8 s over 1 → 8 threads, then is flat. **Python smoothing** keeps creeping down (142.7 → 126.0 s at 16 → 32). Rust is faster at every thread count, **10.4×** at 16 threads.
- **Peak memory:** HKS 9.4–10.8 GiB, KMC 21.7–36.1 GiB; **3.47×** less for HKS at 16 threads (10.07 vs 34.90).
- **Smoothing time and peak memory are identical cold vs warm.** Cache only affects assignment.

## 4. Results: COLD (not plotted, but should be disclosed)

| Threads | HKS assign (Mb/s) | KMC assign (Mb/s) | Extra assignment time vs warm, median (HKS / KMC) |
|---|---|---|---|
| 1 | 6.17 ± 0.18 | 2.19 ± 0.08 | 17 / 56 s |
| 2 | 10.09 ± 0.33 | 4.01 ± 0.10 | 54 / 66 s |
| 4 | 15.96 ± 0.97 | 6.50 ± 0.18 | 41 / 53 s |
| 8 | 24.35 ± 3.51 | 9.67 ± 1.28 | 43 / 62 s |
| 16 | 36.71 ± 5.84 | 11.84 ± 0.30 | 29 / 64 s |
| 32 | 36.67 ± 2.02 | 12.14 ± 1.50 | 42 / 46 s |

- **First-run cost:** a first run adds a roughly fixed **~30–50 s (HKS) / ~45–65 s (KMC)** to assignment, which is the time to read the index from storage. It does not depend on thread count.
- **Effect on scaling:** that fixed cost hides the extra speed at high thread counts. **Cold HKS does not improve past 16 threads** (36.7 vs 36.7 Mb/s), and the HKS/KMC ratio at 16 threads is 3.10× instead of 3.42×.
- **End-to-end at 16 threads, cold:** HKS 1.66 min, KMC 6.61 min (3.98×).
- **Why cold SDs are larger at 8–16 threads:**
  - The cluster's network storage keeps its own cache, which a job cannot clear. The fastest cold HKS runs each started within about 4 minutes of another HKS run, so that cache probably still held the index.
  - The run order was random, so this adds noise but no bias.
  - This uncontrollable noise is one reason the warm numbers are the cleaner measure of compute scaling.

## 5. Advice: how to describe and present this

1. **State the cache condition plainly.** A suggested framing (adapt to the manuscript's voice): each configuration was run twice in succession after clearing the operating-system file cache, and the second run is reported. That run measures the computation with the database already in memory, as for every run after the first on a node. Avoid the jargon "warm" unless it is defined.
2. **Disclose the first-run cost with numbers, in one sentence.** For example: the first run on a node additionally reads the database from storage, adding ~30–50 s (HKS) and ~45–65 s (KMC) to assignment regardless of thread count. Smoothing and peak memory are unaffected. Also note that with that cost included, HKS assignment no longer improves past 16 threads. A reader who runs KaryoScope once should not be surprised by their timing. This keeps the warm-only figure honest.
3. **Delete or replace the old variation paragraph** (main.tex ~line 678: "Run-to-run variation… 129.9 s… first replicate slower than the third… all three replicates are retained"). It describes the confounded old sweep. Point 2 replaces it.
4. **Replicates:** all three replicates are kept at every point and none was dropped; error bars are SD of n = 3.
5. **Correct the claim "Neither step improved past 16 threads."** With warm data it holds for KMC assignment and for Rust smoothing, but not for HKS assignment (+46% from 16 to 32 threads) or, slightly, Python smoothing (−12%).
6. **Table S1 was RE-RUN under the same protocol (job 38434006, done 2026-10-13 22:00). Results and claim-by-claim updates: `$KS_WORK/ks_slurm/tableS1_cold/HANDOFF-tableS1.md`.** The original (job 35247555) ran its 18 runs back-to-back with no cache control, and replicate 1 was the slowest run for every input × backend (8–20%), so its cache state was mixed. Its per-run logs were purged. The re-run reports every Table S1 column for both cold and warm runs. Settled already:
   - **Node:** the original ran on **g-c-1-7-24**. g-c-1-4-02 was a cancelled 64G attempt (job 35247205).
   - **Input sizes (`.fai`):** CHM13 v2.0 3,117,292,070 bp; HG002 v1.1 hap1 **2,947,895,164 bp (2.95 Gb, so Table S1's "2.89 Gbp" is wrong)**; HG002 v1.1 diploid 5,999,424,718 bp.
7. **Same-hardware claim:** "the same hardware as Table S1" is still true (same Ice Lake 2.80 GHz node class).

## 6. Claim-by-claim updates (current text → warm value)

| Location (origin/main `manuscript/main.tex`) | Current | New (warm) |
|---|---|---|
| Caption ~847 | "against the KS_human_CHM13_v2 database" | KS_human_CHM13_v2 (KMC) **and HKS_human_CHM13_v2 (HKS)** |
| Caption 848 | 3.2× at 16 threads (49.1 vs 15.4 Mb/s) | **3.4× (54.8 vs 16.0 Mb/s)** |
| Caption 848 | "Neither backend improved past 16 threads." | KMC did not improve past 16; HKS rose to 79.8 Mb/s at 32 |
| Caption 848 | "The wide interval on HKS at 32 threads comes from one slow replicate." | **delete** |
| Caption 849 | Rust 45.3 → 16.7 s over 1–16 threads, flattening beyond 8 | **42.2 → 13.8 s**, flat from 8 threads (13.8 s at 8, 16 and 32) |
| Caption 849 | 9.2× at 16 (16.7 vs 154.0 s) | **10.4× (13.8 vs 142.7 s)** |
| Caption 850 | HKS 9.4–10.9 GiB; KMC 21.6–36.1 GiB; 3.5× at 16 (10.1 vs 34.9) | HKS **9.4–10.8**; KMC **21.7–36.1**; **3.5× (10.1 vs 34.9)**, effectively unchanged |
| Supp Note 2, ~676 | 15.4 → 49.1 Mb/s at 16 threads (3.2×) | **16.0 → 54.8 (3.4×)** |
| Supp Note 2, ~676 | smoothing 154.0 → 16.7 s at 16 threads (9.2×) | **142.7 → 13.8 s (10.4×)** |
| Supp Note 2, ~676 | "Neither step improved past 16 threads." | see §5 item 5 |
| Supp Note 2, ~676 | 10.1 vs 34.9 GiB at 16 (3.5×) | unchanged (10.1 vs 34.9, 3.5×) |
| Supp Note 2, ~676 | repeat set end-to-end 5.78 → 1.28 min at 16 threads (4.5×) | **5.53 → 1.19 min (4.7×)** |
| Supp Note 2, ~678 | the replicate-variation paragraph | replace per §5 items 2–3 |
| Main text ~196 and Table S1/S2 | 13.7 / 7.9 min etc. | **changed: see `tableS1_cold/HANDOFF-tableS1.md` §6** |

`figures/fig3/data/plot_figS6.py` docstring: it says the data are "extracted
from the existing thread-sweep logs (jobs 35376001 / 35376043) … No new compute was
required". Replace this with job 38373939, the cache protocol above, and the input
file `figS6_phase_timings.warm.tsv`. The plot code itself needs no changes; it was
run unmodified (only `--tsv` / `--out` were given).

## 7. File locations

All under `$KS_WORK/ks_slurm/threadsweep_cold/`:

| File | What it is |
|---|---|
| `run_sweep.sh` | sbatch script, with the design rationale in its header |
| `manifest.tsv` | the 36 points in their random run order |
| `evict.py`, `resident.sh` | cache eviction (`posix_fadvise`) and residency check (`fincore`) |
| `runs/index.tsv` | one row per run (72): backend, threads, rep, mode, node, job, bytes resident before the run, start time, load average, log path |
| `runs/{backend}_t{N}_r{R}_{cold,warm}.log` | full `karyoscope -v` and `/usr/bin/time -v` log for each run |
| `coldsweep_38373939.log` | job log, with the per-run co-tenant count (`others=`) |
| `extract_phases.sh {cold,warm}` | log-to-TSV extractor; reproduces the committed figure rows exactly when run on the old logs |
| **`figS6_phase_timings.warm.tsv`** | **per-run warm data. This replaces `figures/fig3/data/figS6_phase_timings.tsv` in the paper repo.** |
| `figS6_phase_timings.cold.tsv` | per-run cold data |
| **`figs/figS6_warm.{svg,pdf,png}`** | **the warm figure. This replaces `figures/fig3/repeat_speedup_icelake.{svg,pdf,png}`, which main.tex includes as `fig3/repeat_speedup_icelake.pdf`.** |
| `figs/figS6_summary_warm.tsv` | means/SDs. This replaces `figures/fig3/data/figS6_summary.tsv`. |
| `figs/figS6_cold.*`, `figs/figS6_summary_cold.tsv` | the cold versions, for reference |
| `figs/plot_figS6.origin_main.py` | the exact plot script used (copied from paper origin/main) |
| `STATUS.md` | working notes and history |
| `ks-v2.3.0/`, `hks-0.4.0`, `get_featureIDs` | pinned code |

The old sweep (superseded, kept for the record) is in
`$KS_WORK/ks_slurm/threadsweep/`.

To regenerate the figure from the paper repo:
`python figures/fig3/data/plot_figS6.py --tsv <figS6_phase_timings.warm.tsv> --out figures/fig3/repeat_speedup_icelake.svg --summary-tsv figures/fig3/data/figS6_summary.tsv`
