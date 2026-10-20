#!/usr/bin/env python3
"""Pin the ARRAY-COUNT plateau edges exactly, and dump the full qualifying-size list.

The main refinement resolved the repeat_count_main plateau to [1201, 1500].  But
the array count (and the canonical/degraded split) is flatter: it holds 629
through T=1600 and first moves at T=1650.  Because the swept grid was step 50
there, the exact edge is unmeasured -- it is either 1600 or somewhere in
(1600,1650).  Quantization says 1600, but measure it rather than infer it.
Also probes downward for the array-count lower edge, and dumps all qualifying
cluster sizes so the 100-bp quantization claim is on the record in full.
"""
import os as _ks_os
KS_PROJECTS = _ks_os.environ.get("KS_PROJECTS", "$KS_PROJECTS")  # see the repository README, Paths
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import os, sys, shutil, importlib, json, logging
from pathlib import Path

W = Path((KS_WORK + "/a6_sweeps/d4z4_bisect"))
DATA = (KS_PROJECTS + "/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed")
sys.path.insert(0, str(W / "code"))
CONST = "MIN_ARRAY_ELEMENT_SIZE"
_recorded = []


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
    return dict(n_hap=len(rows), arrays=s("num_arrays"), repeat_main=s("repeat_count_main"),
                repeat_total=s("repeat_count_total"), canonical=s("num_canonical"),
                degraded=s("num_degraded")), {r["haplotype"]: int(float(r["num_arrays"])) for r in rows}


def run_point(tag, value, record=False):
    base = W / "points" / tag
    if base.exists():
        shutil.rmtree(base)
    (base / "results").mkdir(parents=True)
    if not (base / "data").exists():
        (base / "data").symlink_to(DATA)
    for m in [m for m in list(sys.modules) if m == "analyze_d4z4_arrays"]:
        del sys.modules[m]
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
    mod = importlib.import_module("analyze_d4z4_arrays")
    setattr(mod, CONST, value)
    if record:
        _recorded.clear()
        original = mod.identify_arrays
        def rec(clusters, gene_regions):
            for st, en in clusters or []:
                fr = mod.get_cluster_annotation_fractions(st, en, gene_regions)
                _recorded.append((en - st, fr["exon"]))
            return original(clusters, gene_regions)
        mod.identify_arrays = rec
    os.environ["D4Z4_BASE_DIR"] = str(base)
    saved = sys.stdout
    try:
        with open(base / "stdout.log", "w") as fh:
            sys.stdout = fh
            mod.main()
    finally:
        sys.stdout = saved
    return read_results(base / "results")


print("=== baseline 1500 (recording all cluster sizes) ===", flush=True)
base_agg, base_hap = run_point("confirm_1500", 1500, record=True)
mod = importlib.import_module("analyze_d4z4_arrays")
thr = mod.EXON_DENSE_THRESHOLD
qual = sorted({sz for sz, ex in _recorded if ex >= thr})
print(f"  baseline {base_agg}", flush=True)
print(f"  all {len(qual)} qualifying sizes: {qual}", flush=True)
non100 = [s for s in qual if s % 100 != 0]
print(f"  qualifying sizes NOT a multiple of 100: {non100 if non100 else 'NONE'}", flush=True)
(W / "qualifying_sizes_full.json").write_text(json.dumps(
    {"exon_dense_threshold": thr, "all_qualifying_sizes": qual,
     "not_multiple_of_100": non100}, indent=2))

rows = []
print("\n=== probe array-count plateau edges ===", flush=True)
for v in [1601, 1700, 1000, 900, 800, 700, 600]:
    agg, hap = run_point(f"confirm_{v}", v)
    ch = sum(1 for h, n in hap.items() if base_hap.get(h) != n)
    d = [k for k in base_agg if base_agg[k] != agg[k]]
    arrays_same = agg["arrays"] == base_agg["arrays"] and \
                  agg["canonical"] == base_agg["canonical"] and \
                  agg["degraded"] == base_agg["degraded"]
    rows.append(dict(value=v, **agg, haps_changed=ch,
                     arrays_and_class_same=int(arrays_same),
                     differs=",".join(d) if d else "none"))
    print(f"  T={v:<6} arrays={agg['arrays']:>5} repeat_main={agg['repeat_main']:>6} "
          f"repeat_total={agg['repeat_total']:>6} canon={agg['canonical']:>4} "
          f"degr={agg['degraded']:>3} haps_changed={ch:>4}  "
          f"array+class {'SAME' if arrays_same else 'CHANGED'}  "
          f"[differs: {','.join(d) if d else 'none'}]", flush=True)

cols = ["value", "n_hap", "arrays", "repeat_main", "repeat_total", "canonical",
        "degraded", "haps_changed", "arrays_and_class_same", "differs"]
(W / "array_count_edges.tsv").write_text("\t".join(cols) + "\n" +
    "\n".join("\t".join(str(r[c]) for c in cols) for r in rows) + "\n")
print("\nwrote array_count_edges.tsv", flush=True)
