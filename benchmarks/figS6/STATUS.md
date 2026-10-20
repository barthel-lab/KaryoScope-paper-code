# Fig. S6 cold-start re-sweep — status / pick-up notes

**Job 38373939** (QOS clong, 48 h limit) STARTED 2026-10-12 15:54 on g-c-1-6-21.
Expected runtime ~14 h (cold ~7.6 h + warm reruns ~6.7 h, estimated from the old
sweep) -> **done ~06:00-07:00 Tue 2026-10-13**. The full cold+warm design is
running (the "warm only for cheap points" trim was discussed but NOT applied).

## Why this sweep exists
User (2026-10-12): every point must be run as a user would see it if it were the
only thing they were running; cache must not be a possible confounder. The old
sweep (jobs 35376001/35376043, `../threadsweep/`) was confounded:
descending thread order put the single cold-index run on t=32 (the 104 s
replicate), rep 3 came from a second job and was 5-30 % faster everywhere, and
both jobs shared their node (scattered CPU affinity, 25 GB of 515 GB free).
`$KS_WORK/figS6_variants/REPLICATE-POLICY.md` (drop-r1 vs keep-all)
is therefore MOOT.

## Design (details in run_sweep.sh header)
- one job, one node, 42 CPUs + 490 GB (cluster rejects --exclusive/--mem=0;
  real submissions cap at -c 42) -> at most 4 CPUs left for other tenants;
  co-tenant count logged per run (`others=` in the job log)
- 36 points = {kmc,hks} x t{1,2,4,8,16,32} x rep{1,2,3}, fixed random order
  (manifest.tsv, seed 20261012)
- before each COLD run: DB + input + binaries + source evicted (evict.py,
  posix_fadvise DONTNEED); fincore residency must be <64 MiB or the job aborts.
  Then an identical WARM rerun.
- pinned: KaryoScope v2.3.0 worktree `ks-v2.3.0/` (PYTHONPATH) + `hks-0.4.0`
  copy + `get_featureIDs` copy (unchanged since v2.3.0)
- outputs to node-local /tmp, deleted between runs
- resumable: rerunning `sbatch run_sweep.sh` skips points already in runs/index.tsv

## Tomorrow: checklist
1. `squeue -u $USER -n ks_coldsweep`; `tail coldsweep_38373939.log` -> "ALL DONE"?
   `wc -l runs/index.tsv` should be 73 (header + 72).
   If it died: just resubmit `run_sweep.sh` (resumes).
2. Sanity: every cold row has resident_bytes_before ~0; `grep others= coldsweep_*.log`
   should be 0 throughout (if not, note which points had co-tenants).
3. `./extract_phases.sh cold > figS6_phase_timings.cold.tsv` (and `warm`).
   Extractor verified to reproduce the committed figure rows exactly on old logs.
4. Plot with KaryoScope-paper `figures/fig3/data/plot_figS6.py` (reads that column
   layout). The paper checkout is SHARED -> use `git worktree`.
5. DECISION FOR USER: plot cold only (what they asked for) or cold + warm.
   Cold includes the index read (~60-90 s HKS), so high-thread HKS throughput
   will drop well below the published 57-71 Mb/s and flatten.
6. Recompute caption numbers in `manuscript/main.tex` (origin/main ~line 847:
   "3.2x at 16 threads (49.1 vs 15.4)", "Neither backend improved past 16",
   "wide interval on HKS at 32 threads comes from one slow replicate" -> delete,
   smoothing 45.3->16.7 s, 9.2x, memory 9.4-10.9 / 21.6-36.1 / 3.5x), and
   Supp Note mentions of the same claims. Also replace the "No new compute" text
   in plot_figS6.py's docstring (it names jobs 35376001/43).

## First data point (sanity, 2026-10-12)
hks t=2 rep2: cold 5:26 wall, warm 4:37 (cold penalty ~49 s). resident 0 before
cold; 10.6 GB resident before warm; loadavg ~2; others=0.

## Early comparison vs the old figure data (2026-10-12, 3.5 points in)
Assignment s, new cold / new warm / old r1,r2,r3:
HKS t2 281.8 / 230.0 / 370.0,290.2,231.1 ; HKS t16 82.0 / 53.5 / 61.3,62.2,57.0 ;
HKS t32 85.7 / 36.7 / 104.0,42.7,40.5 ; KMC t16 256.2 / (pending) / 207.3,198.6,173.0.
=> the old figure is effectively WARM except the 104 s r1 (a cold start, plus extra
slowness probably from the shared node). New warm is 5-15 % faster than old (clean node). Cold adds a
~fixed 30-50 s to HKS lookup (index read). Smoothing unaffected by cache (~14 s).
Cold-only HKS barely gains 16->32 threads (82->86 s); warm still scales (54->37 s).
Told the user already; decision on cold-only vs cold+warm still open.

## RESULT (2026-10-13 10:40) — checklist steps 1-4 DONE
- Job COMPLETED 04:12 (12.3 h); 72/72 rows, exit 0 everywhere, every cold run resident=0, others=0 on all 72.
- `figS6_phase_timings.{cold,warm}.tsv` extracted; figures + summaries in `figs/` (rendered with
  origin/main plot_figS6.py copied to figs/plot_figS6.origin_main.py; paper repo untouched).
- WARM: 32-thread HKS bars gone (SD 0.5 Mb/s); HKS 6.5->54.8->79.8 Mb/s; KMC flat past 16 (16.0/15.5).
  HKS/KMC @16 = 3.42x; smoothing 42.2->13.8 s (1->16), 10.4x; RSS HKS 9.4-10.8, KMC 21.7-36.1, 3.46x.
- COLD: HKS 6.2->36.7 @16 and 36.7 @32 (flat); KMC 12.1 @32. HKS/KMC @16 = 3.10x. Wider bars at t8/t16
  (HKS SD 3.5/5.8 Mb/s): fastest cold HKS runs started <=4 min after another HKS run -> Isilon
  server-side cache (uncontrollable); random order makes it noise, not bias.
- Still open: step 5 (user picks cold / warm / both) and step 6 (caption numbers).

## 2026-10-13: user chose WARM ONLY for the manuscript
Hand-off for the manuscript Claude: `HANDOFF-figS6.md` (results, presentation advice, claim-by-claim numbers, file map). Paper repo still untouched.
