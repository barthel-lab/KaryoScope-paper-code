# A1: comparator runtimes

Wall-clocks for the tools whose runtimes are compared with KaryoScope's, all on
one Intel Xeon Gold 6342 node (Ice Lake, 2.80 GHz) to match Table S1.

| tool | configuration | wall-clock |
|---|---|---|
| RepeatMasker 4.1.6, Dfam 3.8 | `-species human -pa 8` (8 jobs x 4 threads) | **43 h 54 m 35 s**, 97.4 core-hours, 78.9 GiB peak |
| HumAS-HMMER for AnVIL, AS-HORs 3.3.2-120124 | per contig, 16 threads | **53.10 h** whole assembly / **51.76 h** centromeric only, 23/23 contigs both |
| Tiberius, unmasked mammalian model | 1x A100 + 8 CPU | **2 h 15 m 56 s**, 21,981 gene models |
| KaryoScope v2.3.0, all six feature sets | 16 threads | 11.5 min (KMC) / 6.3 min (HKS), *n* = 3 means (Supplementary Table 1) |

## Three things that must be stated with these numbers

1. **The RepeatMasker version is not the one behind KaryoScope's `repeat` feature
   set**, which derives from 4.1.2p1 with the 2022Apr14 Dfam release. These are
   **runtimes only**, with no concordance computed.
2. **`-pa` parallelises the batch search only.** `ProcessRepeats` is
   single-threaded and dominates the wall clock (average CPU use was 2.2 of 32
   allocated cores). Do not read 43.9 h as a parallel-efficiency figure.
3. **Extracting the centromeres does not reduce HumAS-HMMER's runtime.** It took
   53.10 h over 2,745 Mb compared with 51.76 h over 271 Mb (10x less sequence,
   3% less time), because nhmmer's cost tracks alpha-satellite *content*, not
   total bases. HOR (higher-order repeat) call counts are identical between the
   two modes, so the extraction itself is correct.

## The nhmmer alphabet trap

Contigs chr13, chr21, and chr22 fail with "Unable to guess alphabet" in 0.2 s.
Those three HG002 hap1 contigs open with a pure telomeric `(TAACCC)n` tract with
**no G** for longer than Easel's ~4 kb guess window (first G at residue 16,721 /
9,818 / 9,433). With only A, C, and T visible, the guesser cannot separate DNA
from protein.

**Fix: `nhmmer --dna`** (asserted, so no cost change). `--tformat fasta` does not
help. The fix is patched into the vendored `HumAS-HMMER_for_AnVIL/hmmer-run.sh`.

The failure was invisible because `run_a1b_array.sh` captured `$RC` but never
exited on it, so SLURM logged COMPLETED. `exit $RC` added.

## Not done, deliberately

`RM2Bed.py` needs pandas, which the RepeatMasker environment lacks. The `.out`
file is intact, and the conversion takes minutes offline. Because of the version
caveat above, the converted file would not be used for concordance, so the
conversion was not run.
