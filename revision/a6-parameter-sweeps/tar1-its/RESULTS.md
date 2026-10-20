# A6 §6 — TAR1 / ITS merge distance: results

**Run:** 2026-08-24, Intel Xeon Gold 6342 (Ice Lake) @ 2.80 GHz, the Table S1 node class.
**Jobs:** `36897594` (inputs), `36900174` (sweep).
**Verdict:** **sweep-backed row — but the evidence is a bridging ceiling, not the knee the brief expected.**

---

## 1. What was swept

`main.tex:338` (TAR1) and `:340` (ITS) both use a 1 kb merge distance with one
stated rationale: consolidate RepeatMasker's fragmented per-element hits into one
feature per array.

Source is the S3 RepeatMasker release the subtelomeric set was built from —
`chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed`, 342,862,439 bytes,
sha256 `abd319705cda7090d0c5b6ac002f0af6538ef22217965b38135f9185355fd285`.
**This is a different file from the UCSC bigBed the `repeat` set uses**; the recipe
is explicit that the two representations are not interchangeable.

Column layout, verified rather than assumed:
`chrom start end name score strand class family divergence id`

* **TAR1** — field 4 contains `TAR1`. 201 hits, 194 after collapsing exact overlaps.
* **ITS** — field 7 is `Simple_repeat` **and** field 4 contains a rotation of
  TTAGGG or CCCTAA. 337 candidates, then those overlapping a `seqtk telo` region
  are dropped, leaving 289. `seqtk telo` on CHM13 v2 gives 48 regions, 132,815 bp.

## 2. The sweep

| feature | merge_d | intervals | total bp | median | max |
|---|---|---|---|---|---|
| TAR1 | 0 | 194 | 135,875 | 575 | 2,742 |
| TAR1 | 100 | 185 | 136,187 | 618 | 2,742 |
| TAR1 | 500 | 181 | 137,521 | 629 | 2,742 |
| **TAR1** | **1000** | **180** | **138,396** | **625** | **2,742** |
| TAR1 | 2000 | 175 | 143,448 | 629 | 2,742 |
| TAR1 | 5000 | 175 | 143,448 | 629 | 2,742 |
| TAR1 | 10000 | 174 | 150,598 | 629 | **9,131** |
| ITS | 0 | 286 | 52,723 | 52 | 5,362 |
| ITS | 100 | 275 | 53,055 | 52 | 5,362 |
| ITS | 500 | 270 | 54,364 | 55 | 5,362 |
| **ITS** | **1000** | **262** | **59,889** | **55** | **5,362** |
| ITS | 2000 | 260 | 63,010 | 55 | 5,362 |
| ITS | 5000 | 255 | 84,647 | 56 | **11,154** |
| ITS | 10000 | 254 | 90,569 | 55 | 11,154 |

### There is no knee

The brief predicted fragment count falling steeply then flattening, with 1 kb past
the knee. **That is not what the data does.** Across four orders of magnitude the
interval count falls only 194 → 174 for TAR1 (10%) and 286 → 254 for ITS (11%), and
the decline is shallow throughout rather than steep-then-flat. RepeatMasker's TAR1
and ITS hits are simply not heavily fragmented at the sub-kilobase scale, so there
is little for a consolidation threshold to consolidate.

This is the same shape as the `max_gap` result: the parameter is far less
consequential than the manuscript's framing implies.

### What the threshold *does* control is bridging

Total bases labelled is the informative column, and it is not flat. Merging spans
the gaps it bridges, so total bases can only grow:

| | at d=1000 vs d=0 | at d=10000 vs d=0 |
|---|---|---|
| TAR1 | +2,521 bp (**+1.9%**) | +14,723 bp (+10.8%) |
| ITS | +7,166 bp (**+13.6%**) | +37,846 bp (+71.8%) |

At 1 kb the absolute inflation is small — 14 TAR1 merges bridging ~180 bp of gap
each, 24 ITS merges bridging ~300 bp each. Those are plausibly one fragmented
element apiece. The ITS percentage looks alarming only because the ITS set is
small (52.7 kb total).

**The ceiling is visible and it is above 1 kb.** The longest interval is stable at
2,742 bp (TAR1) and 5,362 bp (ITS) from d=0 all the way to d=2000, then jumps to
9,131 bp and 11,154 bp. That is the signature of the merge bridging genuinely
separate elements rather than reassembling one. For ITS it happens by d=5000, for
TAR1 by d=10000.

So 1 kb sits in a wide safe band: above the gaps that fragment single elements,
and comfortably below where distinct elements start fusing.

## 3. The filter/merge ordering does not matter here

The brief flagged that the `seqtk telo` filter must be applied before merging, so
an interior hit cannot bridge into a terminal array. Tested directly at d=1000:

```
filter-then-merge: 262 intervals, 59,889 bp
merge-then-filter: 262 intervals, 59,889 bp
```

**Identical.** The 48 telomere-overlapping candidates are far enough from interior
hits that a 1 kb merge never bridges them. The ordering is still the correct one to
write down and the code should keep it, but it is not load-bearing at the published
threshold. It could become so at larger *d*; that was not tested beyond 1 kb.

## 4. Recommended row

**Verdict: sweep-backed**, with the justification restated. The right sentence is
not "1 kb sits past the knee in fragment count" — there is no knee. It is:

> Merge distance is set to 1 kb. Interval count is insensitive to it (TAR1
> 194 -> 174 and ITS 286 -> 254 across 0-10 kb), so the threshold is chosen by
> what it absorbs rather than what it consolidates: at 1 kb it adds 1.9% (TAR1)
> and 13.6% (ITS) to the total labelled span, bridging gaps of a few hundred bases
> within fragmented elements, while the longest interval is unchanged from the
> unmerged set. Above 2 kb the longest interval more than doubles, indicating that
> distinct elements begin to fuse.

That is a stronger claim than the current one because it names the failure mode the
threshold is avoiding and shows where it starts.

## 5. Files

```
a6_sweeps/tar1_its/run_fetch.sh, run_sweep.sh
a6_sweeps/tar1_its/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed   (327 MB)
a6_sweeps/tar1_its/chm13v2.0.telo.bed
a6_sweeps/tar1_its/results/tar1_its_sweep.tsv                     <- the table above
a6_sweeps/tar1_its/results/{tar1,its}.raw.bed, {tar1,its}.merged_d*.bed
```

## 6. Not done

Ordering was checked only at d=1000. The per-interval size distribution is in the
TSV (min/q1/median/q3/max) but was not plotted.
