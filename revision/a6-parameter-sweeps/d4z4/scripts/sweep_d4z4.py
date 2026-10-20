#!/usr/bin/env python3
"""A6 Sweep 4 (D4Z4) -- the merge thresholds in HPRC-D4Z4:analyze_d4z4_arrays.py.

Every swept constant is a MODULE-LEVEL global read at call time, so the driver
overrides the module attribute and calls main() -- Floris's logic is untouched.
The only source edit is `base_dir`, made env-overridable in one line because the
upstream value is a macOS path.

Varies one constant at a time around its published value, over chr4 + chr10
(the only two the script processes), 582 haplotypes. ~16 s per point.

Reported: total arrays and total repeat units -- the repeat count is the
published quantity, so that is the number whose stability matters. Also
canonical/degraded split and the number of haplotypes whose array count changes,
because an aggregate can hide compensating moves.

Groups:
  7a  the seven manuscript-stated thresholds
  7b  EXON_DENSE_THRESHOLD -- undisclosed, but it defines the primitive
      (arrays are runs of exon-dense BSR clusters) so it sits upstream of 7a
  7c  BSR_MERGE_THRESHOLD, MIN_ARRAY_ELEMENT_SIZE -- undisclosed; swept to
      settle whether they belong in the table
"""
import os as _ks_os
KS_PROJECTS = _ks_os.environ.get("KS_PROJECTS", "$KS_PROJECTS")  # see the repository README, Paths
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import os, sys, shutil, importlib, json
from pathlib import Path

W = Path((KS_WORK + "/a6_sweeps/d4z4"))
DATA = (KS_PROJECTS + "/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed")
sys.path.insert(0, str(W / "code"))

# constant -> (published, [values], group)
SWEEP = {
    "ARRAY_MERGE_GAP":          (4000,  [1000, 2000, 3000, 4000, 6000, 8000],       "7a"),
    "NOVEL_MERGE_THRESHOLD":    (0.7,   [0.5, 0.6, 0.7, 0.8, 0.9],                  "7a"),
    "EXON_MERGE_THRESHOLD":     (0.6,   [0.4, 0.5, 0.6, 0.7, 0.8],                  "7a"),
    "COMBINED_MERGE_THRESHOLD": (0.8,   [0.6, 0.7, 0.8, 0.9, 1.0],                  "7a"),
    "BOUNDARY_WINDOW_SIZE":     (3000,  [1000, 2000, 3000, 4000, 5000],             "7a"),
    "BOUNDARY_BLOCK_THRESHOLD": (0.3,   [0.1, 0.2, 0.3, 0.4, 0.5],                  "7a"),
    "ARRAY_GAP_THRESHOLD":      (10000, [5000, 7500, 10000, 15000, 20000],          "7a"),
    "EXON_DENSE_THRESHOLD":     (0.5,   [0.3, 0.4, 0.5, 0.6, 0.7],                  "7b"),
    "BSR_MERGE_THRESHOLD":      (400,   [200, 300, 400, 600, 800],                  "7c"),
    "MIN_ARRAY_ELEMENT_SIZE":   (1500,  [1000, 1250, 1500, 1750, 2000],             "7c"),
    # the remaining undisclosed constants from the inventory, swept so the table
    # can say what was measured rather than what was not looked at
    "PLAM_MAX_EXON_FRACTION":   (0.35,  [0.25, 0.30, 0.35, 0.40, 0.45],             "7c"),
    "WINDOW_SIZE":              (250000,[150000, 200000, 250000, 350000, 500000],   "7c"),
    "TRUNCATION_THRESHOLD":     (1000,  [500, 750, 1000, 1500, 2000],               "7c"),
    "TRACT_GAP_THRESHOLD":      (10000, [5000, 7500, 10000, 15000, 20000],          "7c"),
}


def read_results(results_dir):
    f = Path(results_dir) / "d4z4_array_analysis.tsv"
    rows, hdr = [], None
    for i, line in enumerate(f.read_text().splitlines()):
        parts = line.split("\t")
        if i == 0:
            hdr = parts
            continue
        rows.append(dict(zip(hdr, parts)))
    def s(col):
        return sum(int(float(r[col])) for r in rows if r.get(col) not in (None, "", "nan"))
    per_hap = {r["haplotype"]: int(float(r["num_arrays"])) for r in rows}
    return dict(n_hap=len(rows), arrays=s("num_arrays"),
                repeat_main=s("repeat_count_main"), repeat_total=s("repeat_count_total"),
                canonical=s("num_canonical"), degraded=s("num_degraded")), per_hap


def run_point(tag, overrides):
    base = W / "points" / tag
    if base.exists():
        shutil.rmtree(base)
    (base / "results").mkdir(parents=True)
    link = base / "data"
    if not link.exists():
        link.symlink_to(DATA)

    # Fresh import per point so no state leaks between runs.
    for m in [m for m in list(sys.modules) if m == "analyze_d4z4_arrays"]:
        del sys.modules[m]
    import logging
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
    mod = importlib.import_module("analyze_d4z4_arrays")
    for k, v in overrides.items():
        assert hasattr(mod, k), f"{k} is not a module attribute -- monkey-patch would silently do nothing"
        setattr(mod, k, v)
    os.environ["D4Z4_BASE_DIR"] = str(base)
    out = base / "stdout.log"
    saved = sys.stdout
    try:
        with open(out, "w") as fh:
            sys.stdout = fh
            mod.main()
    finally:
        sys.stdout = saved
    return read_results(base / "results")


def main():
    rows = []
    print("=== baseline (all published values) ===", flush=True)
    base_agg, base_per_hap = run_point("sweep_baseline", {})
    print(f"  {base_agg}", flush=True)

    for const, (published, values, group) in SWEEP.items():
        assert published in values, f"{const}: published {published} not in {values}"
        for v in values:
            tag = f"{const}_{v}"
            agg, per_hap = run_point(tag, {const: v})
            changed = sum(1 for h, n in per_hap.items() if base_per_hap.get(h) != n)
            rows.append(dict(group=group, constant=const, value=v,
                             is_published=int(v == published), **agg,
                             haps_changed=changed,
                             d_arrays=agg["arrays"] - base_agg["arrays"],
                             d_repeat_main=agg["repeat_main"] - base_agg["repeat_main"],
                             d_repeat_total=agg["repeat_total"] - base_agg["repeat_total"]))
            print(f"  {tag:<38} arrays={agg['arrays']:>5} repeat_main={agg['repeat_main']:>6} "
                  f"repeat_total={agg['repeat_total']:>6} haps_changed={changed:>4}"
                  + ("  <- published" if v == published else ""), flush=True)

    cols = ["group", "constant", "value", "is_published", "n_hap", "arrays",
            "repeat_main", "repeat_total", "canonical", "degraded",
            "haps_changed", "d_arrays", "d_repeat_main", "d_repeat_total"]
    dest = W / "d4z4_param_sweep.tsv"
    dest.write_text("\t".join(cols) + "\n" +
                    "\n".join("\t".join(str(r[c]) for c in cols) for r in rows) + "\n")
    (W / "d4z4_baseline.json").write_text(json.dumps(base_agg, indent=2))
    print(f"\nwrote {dest} ({len(rows)} rows)")


main()
