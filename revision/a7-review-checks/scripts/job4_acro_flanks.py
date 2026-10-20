#!/usr/bin/env python3
"""Job 4 -- what flanks the long multilabel runs in the acrocentric set, and did smoothing change them.

Run definition: a maximal stretch of consecutive presmoothed intervals whose labels are NOT leaves
(internal nodes of the acrocentric hierarchy, or `novel`), bounded on each side by a leaf-labelled
interval (or the sequence end). Because annotations tile, the flanking leaf interval is directly
adjacent to the run. A run is reported when its internal-node (multilabel) bp exceeds 1 kb; novel
bp inside the run is counted separately. run_label = the internal label with the most bp in the run.
Expected outcome from the flank rule: reassign iff LCA(left flank, right flank) is a strict
descendant of run_label. Observed = smoothed labels over the run's internal-labelled bases."""
import os as _ks_os
KS_SOFTWARE = _ks_os.environ.get("KS_SOFTWARE", "$KS_SOFTWARE")  # see the repository README, Paths
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import sys, os, bisect, collections
sys.path.insert(0, os.path.dirname(__file__))
from common_mask import *
HIER = (KS_SOFTWARE + "/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv")
A4 = (KS_WORK + "/a4_smoothing")
RES = (KS_WORK + "/a7_review_checks/results")
INPUTS = {"chm13": (f"{A4}/annot_chm13/chm13v2.0.HKS_human_CHM13_v2.acrocentric.presmoothed.bed", f"{A4}/annot_chm13/chm13v2.0.HKS_human_CHM13_v2.acrocentric.smoothed.bed"),
          "hg002hap1": (f"{A4}/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.acrocentric.presmoothed.bed", f"{A4}/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.acrocentric.smoothed.bed")}
MIN_INTERNAL_BP = 1000
MAX_GAP = 1000
nodes, leaves, parent_of = parse_hierarchy(HIER)["acrocentric"]
lvl = make_level(nodes, leaves)
print("acrocentric hierarchy leaves:", sorted(leaves), " internal:", sorted(nodes - leaves))

def by_chrom(path):
    cur = None; buf = []
    for line in open(path):
        c, s, e, lab = line.rstrip("\n").split("\t")[:4]
        if c != cur:
            if cur is not None: yield cur, buf
            cur, buf = c, []
        buf.append((int(s), int(e), lab))
    if cur is not None: yield cur, buf

rows = []; stats = collections.defaultdict(collections.Counter); allruns = collections.Counter()
for asm, (pre_f, post_f) in INPUTS.items():
    assert "presmoothed" in pre_f and "presmoothed" not in os.path.basename(post_f).replace("acrocentric.smoothed", "")
    for (c1, P), (c2, S) in zip(by_chrom(pre_f), by_chrom(post_f)):
        assert c1 == c2, (c1, c2)
        if c1 == "chrM": continue
        s_starts = [x[0] for x in S]
        def smoothed_over(s, e):
            out = collections.Counter(); k = bisect.bisect_right(s_starts, s) - 1
            while k < len(S) and S[k][0] < e:
                lo, hi = max(s, S[k][0]), min(e, S[k][1])
                if hi > lo: out[S[k][2]] += hi - lo
                k += 1
            return out
        n = len(P); i = 0
        while i < n:
            if lvl(P[i][2]) == "leaf": i += 1; continue
            j = i
            while j < n and lvl(P[j][2]) != "leaf": j += 1
            run = P[i:j]
            internal = collections.Counter(); novel_bp = 0
            for s, e, lab in run:
                if lvl(lab) == "internal": internal[lab] += e - s
                else: novel_bp += e - s
            ibp = sum(internal.values())
            allruns["runs_with_internal_bp"] += 1 if ibp else 0
            if ibp > MIN_INTERNAL_BP:
                run_label = max(internal, key=internal.get)
                left = P[i-1] if i > 0 else None; right = P[j] if j < n else None
                ll = left[2] if left else "SEQ_START"; rl = right[2] if right else "SEQ_END"
                if left and right:
                    L = lca(parent_of, ll, rl)
                    expected = "reassign" if (L != run_label and run_label in ancestors(parent_of, L)) else "unchanged"
                else: L = "NA"; expected = "NA"
                # EFFECTIVE flanks, mirroring the window scan: walk outward past intervals labelled with the run
                # label itself, past novel, and past unrelated-branch intervals (cumulative <= MAX_GAP bp), and
                # stop at (a) a leaf inside the run label's subtree, or (b) a label that is a strict ancestor of
                # the run label (caps the window: LCA can be no more specific than it). If the unrelated span
                # exceeds MAX_GAP the window breaks and the boundary is the run label itself.
                def eff_flank(start_idx, step):
                    k = start_idx; skipped_unrel = 0
                    while 0 <= k < n:
                        s_, e_, lab = P[k]
                        if lab.lower() == "novel": k += step; continue
                        anc_of_lab = ancestors(parent_of, lab)
                        if lab == run_label: k += step; continue
                        if run_label in anc_of_lab:            # descendant of run label
                            if lab in leaves: return lab, skipped_unrel
                            k += step; continue                 # a more specific internal node: pass through
                        if lab in ancestors(parent_of, run_label):  # strict ancestor -> caps the window
                            return lab, skipped_unrel
                        skipped_unrel += e_ - s_                # unrelated branch
                        if skipped_unrel > MAX_GAP: return run_label, skipped_unrel
                        k += step
                    return "SEQ_END", skipped_unrel
                effL, skL = eff_flank(i - 1, -1); effR, skR = eff_flank(j, +1)
                if effL in ("SEQ_END",) or effR in ("SEQ_END",): effLCA = "NA"; expected_eff = "NA"
                else:
                    effLCA = lca(parent_of, effL, effR)
                    expected_eff = "reassign" if (effLCA != run_label and run_label in ancestors(parent_of, effLCA)) else "unchanged"
                # observed: smoothed labels over the internal-labelled intervals only
                sm = collections.Counter(); changed = 0
                for s, e, lab in run:
                    if lvl(lab) != "internal": continue
                    d = smoothed_over(s, e)
                    for k2, v in d.items():
                        sm[k2] += v
                        if k2 != lab: changed += v
                observed = "unchanged" if changed == 0 else ("reassigned_fully" if changed == ibp else "partial")
                consistent = "NA" if expected == "NA" else str((expected == "reassign") == (observed != "unchanged"))
                consistent_eff = "NA" if expected_eff == "NA" else str((expected_eff == "reassign") == (observed != "unchanged"))
                stats[asm][(expected, observed)] += 1; stats[asm + "_eff"][(expected_eff, observed)] += 1
                rows.append((asm, c1, run[0][0], run[-1][1], run[-1][1] - run[0][0], ibp, novel_bp, run_label,
                             ";".join(f"{k2}:{v}" for k2, v in internal.most_common()), len(internal) == 1,
                             ll, left[1] - left[0] if left else "", rl, right[1] - right[0] if right else "", L, expected,
                             observed, changed, ";".join(f"{k2}:{v}" for k2, v in sm.most_common()), consistent,
                             effL, skL, effR, skR, effLCA, expected_eff, consistent_eff))
            i = j
        print(f"  {asm} {c1}: done ({len(rows)} rows so far)", flush=True)
hdr = ["assembly", "chromosome", "run_start", "run_end", "span_bp", "multilabel_bp", "novel_bp_in_run", "run_label", "labels_in_run", "single_label_run",
       "left_flank_label", "left_flank_bp", "right_flank_label", "right_flank_bp", "flank_lca", "expected_by_flank_rule", "observed_after_smoothing", "changed_bp", "smoothed_labels_over_run", "rule_consistent", "eff_left_flank", "eff_left_skipped_unrelated_bp", "eff_right_flank", "eff_right_skipped_unrelated_bp", "eff_flank_lca", "expected_by_effective_rule", "effective_rule_consistent"]
with open(f"{RES}/acro_multilabel_flanks.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows: o.write("\t".join(str(x) for x in r) + "\n")
print("\n=== summary: runs with multilabel bp > 1 kb ===")
for asm in INPUTS:
    tot = sum(stats[asm].values()); print(f"{asm}: {tot} runs")
    for (exp, obs), n in sorted(stats[asm].items()): print(f"   expected={exp:9} observed={obs:16} n={n}")
    by_label = collections.Counter(r[7] for r in rows if r[0] == asm); print("   run labels:", dict(by_label))
    incons = [r for r in rows if r[0] == asm and r[19] == "False"]
    print(f"   immediate-flank rule inconsistent runs: {len(incons)}")
    print(f"   EFFECTIVE-flank rule:")
    for (exp, obs), cnt in sorted(stats[asm + "_eff"].items()): print(f"      expected={exp:9} observed={obs:16} n={cnt}")
    incons2 = [r for r in rows if r[0] == asm and r[-1] == "False"]
    print(f"   effective-flank rule inconsistent runs: {len(incons2)}")
    for r in incons2[:40]: print("     ", r[1], r[2], r[3], "label", r[7], "imm flanks", r[10], r[12], "eff flanks", r[20], r[22], "efflca", r[24], "expected", r[25], "observed", r[16], r[18])
