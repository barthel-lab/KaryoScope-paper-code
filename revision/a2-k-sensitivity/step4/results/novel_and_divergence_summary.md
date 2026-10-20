# A2 Step 4, Jobs 1 and 2 -- novel fraction vs k, and the divergence series at k = 21/31/51

Compute 2026-09-03. Ice Lake 2.80 GHz nodes, karyoscope 2.3.0, HKS 0.4.0, 16 threads. Databases
`HKS_human_CHM13_v2_k{11,21,31,41,51,61}` (the six-k family built 2026-08-05 from identical inputs).
Deliverables: `novel_by_k.tsv` (144 rows), `divergence_by_k.tsv` (24 rows). Provenance per run in
`job1_novel/results/prov_*.txt` and `job2_divergence/*/results/provenance.txt`.

## Job 1 -- novel fraction vs k (no re-annotation; the A2 Step 2 BEDs, job 35382788 / 35382789)

Method = A4 novel-v3 (N runs masked out with the k-1 left extension fill_gaps applies, chrM excluded,
a run = coordinate-adjacent `novel` intervals). Anchor: HG002 hap1 k=31 gives 75,396,753 novel bp of
2,940,260,101 evaluated = 2.5643%, 86.90% in runs <= k, 0.60% in runs > 1 kb, all exact to the A4 and
Methods figures. Pre-smoothing novel bp is identical across the six feature sets at every k (it is a
property of the base index), so one table serves all sets; post-smoothing differs per set (in the TSV).

HG002 hap1 (repeat set shown; presmoothed novel identical for all sets):

| k | novel pre (bp) | novel pre % | novel post % | runs <= k | runs <= 31 bp | runs > 1 kb | median run | longest run |
|---|---|---|---|---|---|---|---|---|
| 11 | 11 | 0.0000 | 0.0000 | 100% | 100% | 0 | -- | -- |
| 21 | 42,662,055 | 1.4510 | 0.0540 | 91.07% | 94.93% | 0.008% | 16 | 1,450 |
| 31 | 75,396,753 | 2.5643 | 0.0406 | 86.90% | 86.90% | 0.605% | 31 | 9,673 |
| 41 | 108,712,284 | 3.6974 | 0.0700 | 83.27% | 13.87% | 0.581% | 41 | 13,647 |
| 51 | 141,298,826 | 4.8057 | 0.1037 | 79.84% | 7.15% | 0.654% | 51 | 13,702 |
| 61 | 172,642,557 | 5.8717 | 0.1428 | 76.52% | 3.53% | 0.740% | 61 | 13,742 |

CHM13 (self-annotation): 0 novel bp at every k, pre and post smoothing.

What moved and what did not. The pre-smoothing novel fraction grows almost linearly with k, about
0.11 percentage points per unit of k (1.45% at 21 to 5.87% at 61), as expected when each divergent
site breaks k k-mers; the median novel run is exactly k for k >= 31, confirming that the bulk of
novel sequence is single-site divergence. The share of novel bases in runs <= k falls slowly with k
(91% to 77%). The fraction in runs > 1 kb -- the defensible "candidate new sequence" figure -- is
essentially flat at 0.58-0.74% for k = 31-61, but is only 0.008% at k = 21 (longest run 1.45 kb):
at k = 21 the index is dense enough that genuinely absent sequence still finds 21-mer hits.
Post-smoothing novel is minimal at k = 31 (0.041%) and rises to 0.14% at k = 61; k = 21 is 0.054%.
k = 11 is the degenerate control (99.997% of 11-mer space is in the index): 11 novel bp genome-wide.

## Job 2 -- divergence series at k = 21, 31, 51 (A5 design; mutate_fasta.py seed 42)

The four mutated CHM13 genomes were regenerated and verified byte-identical to A5's on `.fai`,
mutation summary and projected truth. Concordance = the published chain (`ks_concord.sh`,
`--kmer-size K`, repeat set with the tandem_repeat merge; region set also run). Control: the k = 31
rows from `HKS_human_CHM13_v2_k31` are IDENTICAL to A5's rows from the shipped `HKS_human_CHM13_v2`
in all 8 cells (4 conditions x 2 sets), to the last digit.

Repeat set (novel = presmoothed novel fraction of evaluated bases; prediction = 1-(1-r)^k):

| condition | k | novel measured | novel predicted | exact pre | exact post | incompatible post |
|---|---|---|---|---|---|---|
| sub 0.001 | 21 | 0.0171 | 0.0208 | 80.34% | 96.23% | 1.84% |
| sub 0.001 | 31 | 0.0277 | 0.0305 | 87.17% | 98.29% | 0.70% |
| sub 0.001 | 51 | 0.0480 | 0.0497 | 91.29% | 98.98% | 0.27% |
| sub 0.01 | 21 | 0.1579 | 0.1903 | 67.62% | 92.70% | 2.34% |
| sub 0.01 | 31 | 0.2453 | 0.2677 | 66.87% | 95.50% | 0.99% |
| sub 0.01 | 51 | 0.3893 | 0.4010 | 58.29% | 94.46% | 0.44% |
| sub 0.05 | 21 | 0.5740 | 0.6594 | 31.37% | 81.65% | 3.83% |
| sub 0.05 | 31 | 0.7583 | 0.7961 | 20.21% | 83.84% | 1.82% |
| sub 0.05 | 51 | 0.9188 | 0.9269 | 7.58% | 66.41% | 3.60% |
| indel 0.005 | 21 | 0.0907 | 0.0999 | 74.12% | 94.29% | 2.05% |
| indel 0.005 | 31 | 0.1392 | 0.1439 | 77.13% | 96.99% | 0.77% |
| indel 0.005 | 51 | 0.2233 | 0.2256 | 74.48% | 97.14% | 0.30% |

Region set post-smoothing exact agreement (k = 21 / 31 / 51): sub 0.001 87.9 / 93.1 / 97.5%;
sub 0.01 85.6 / 91.9 / 96.6%; sub 0.05 73.7 / 88.0 / 84.9%; indel 0.005 86.9 / 92.8 / 97.3%.

What moved and what did not. (1) The measured novel fraction follows 1-(1-r)^k at every k and rate,
always a little below the prediction (k-mers rescued by other copies of a repeat); the brief's
expected 0.19 / 0.27 / 0.40 at r = 0.01 is measured as 0.158 / 0.245 / 0.389. (2) Pre-smoothing
agreement falls with k at 1% and 5% divergence, as expected. (3) Post-smoothing agreement is much
flatter, and the k-optimum depends on divergence: at 0.1% substitution and at the indel rate,
k = 51 is marginally best; at 1% k = 31 is best (95.5%) with k = 51 1.0 point behind and k = 21
2.8 points behind; at 5% k = 51 collapses (66.4% vs 83.8% at k = 31), because 92% of 51-mers are
broken and smoothing has too little to impute from. So a larger k buys accuracy only while the query
stays within ~1% of the reference; k = 31 is the value that is never worse than 1.7 points off the
best at any tested divergence. (4) The region set's incompatible fraction shrinks monotonically with
k at every divergence <= 1% (11.7% -> 6.6% -> 2.3% at r = 0.001): the arm-versus-centromere floor is
a low-k shared-k-mer effect, not a divergence effect.
