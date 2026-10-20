#!/usr/bin/env python3
"""Job 5 -- same-contig table for the multi-array (and mixed / degraded) D4Z4 haplotypes.
Source: the frozen k=31 caller output that reproduced the census bit-for-bit (handoff #6):
a2_step4/job3_d4z4/base_published/results/d4z4_array_analysis.tsv (== base_v1like_k31, cmp-identical)."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import csv, json, sys
SRC = (KS_WORK + "/a2_step4/job3_d4z4/base_published/results/d4z4_array_analysis.tsv")
WIN = (KS_WORK + "/a2_step4/common/d4z4_windows.tsv")
OUT = (KS_WORK + "/a7_review_checks/results/d4z4_multiarray_contigs.tsv")
win = {}
for line in open(WIN):
    p = line.rstrip("\n").split("\t")
    if p[0] == "name" or len(p) < 8: continue
    win[p[0]] = dict(orig=p[1], chrom=p[2], hap=p[3], fasta=p[4], contig_len=int(p[5]), win_start=int(p[6]), win_end=int(p[7]))
rows = list(csv.DictReader(open(SRC), delimiter="\t"))
assert len(rows) == 582, len(rows)
out = []; n_hap = 0
for r in rows:
    na, nc, nd = int(r["num_arrays"]), int(r["num_canonical"]), int(r["num_degraded"])
    if na < 2 and nd < 1: continue
    n_hap += 1
    cls = "multi-array" if na >= 2 else "single"
    if nc > 0 and nd > 0: cls += "+mixed"
    elif nd > 0 and nc == 0: cls += "+degraded-only"
    hap = r["haplotype"]; w = win[hap]
    contig = hap.split("#")[-1]; orig = w["orig"]
    arrays = json.loads(r["arrays_info"]) if r["arrays_info"] else []
    assert len(arrays) == na, (hap, na, len(arrays))
    arrays.sort(key=lambda a: a["start"])
    contigs = {orig}  # every array of a haplotype comes from the one contig the caller was given
    prev_end = None
    for a in arrays:
        term = "pLAM" if a.get("has_plam") else ("ct_terminal" if a.get("has_ct_terminal") else "none")
        gap = "" if prev_end is None else a["start"] - prev_end
        out.append([r["sample"], r["chromosome"], r["hap_num"], hap, orig, contig, cls, na, nc, nd, r["terminal_feature"], r["notes"],
                    a["index"], a["type"], a["start"], a["end"], w["win_start"] + a["start"], w["win_start"] + a["end"],
                    a["element_count"], term, a.get("plam_size", ""), a.get("plam_end", ""), gap, "TRUE", len(contigs)])
        prev_end = a["end"]
hdr = ["sample", "chromosome", "haplotype", "haplotype_id", "contig_original_name", "contig_accession", "class", "num_arrays", "num_canonical", "num_degraded",
       "terminal_feature_haplotype", "caller_notes", "array_index", "array_type", "array_start_window", "array_end_window", "array_start_contig", "array_end_contig",
       "element_count", "array_terminal_label", "plam_size", "plam_end_window", "gap_from_previous_array_bp", "all_arrays_same_contig", "n_distinct_contigs"]
with open(OUT, "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in out: o.write("\t".join(str(x) for x in r) + "\n")
print(f"haplotypes written: {n_hap}; array rows: {len(out)}")
import collections
print(collections.Counter(r[6] for r in out if r[12] == 1))
