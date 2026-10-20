#!/usr/bin/env python3
"""A6 Sweep 4 refinement -- resolve the MIN_ARRAY_ELEMENT_SIZE plateau EXACTLY.

The step-250 sweep established only that the first change above the published
1500 lies somewhere in (1500, 1750].  That is the grid resolution, not the
margin.  This driver resolves both plateau edges to the base pair.

Method.  MIN_ARRAY_ELEMENT_SIZE enters at exactly one place
(analyze_d4z4_arrays.py:169):

    is_exon_dense = (fractions['exon'] >= EXON_DENSE_THRESHOLD and
                     size >= MIN_ARRAY_ELEMENT_SIZE)

so the set of exon-dense clusters is a step function of the threshold whose
breakpoints are the observed cluster sizes among clusters that pass the exon
test.  Those sizes are data, not a grid.  Phase 1 records them by wrapping
identify_arrays (module-global lookup, same monkey-patch seam the sweep uses --
Floris's logic still runs unmodified).  With Q = the qualifying sizes:

    s0 = min{s in Q : s >= 1500}   -> exon-dense set is identical for T in [1500, s0]
    t0 = max{s in Q : s <  1500}   -> exon-dense set is identical for T in (t0, 1500]

Counts can only move where the set moves, so the reported plateau is a superset
of [t0+1, s0] and its edges lie on those same breakpoints.  Phase 3 walks the
breakpoints outward, running the real pipeline at each, until the reported
counts actually differ -- so the answer needs no monotonicity assumption and no
step size.  Phase 2 also runs a step-50 scan over [1000, 1750] for the record.
"""
import os as _ks_os
KS_PROJECTS = _ks_os.environ.get("KS_PROJECTS", "$KS_PROJECTS")  # see the repository README, Paths
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import os, sys, shutil, importlib, json, logging
from pathlib import Path

W = Path((KS_WORK + "/a6_sweeps/d4z4_bisect"))
DATA = (KS_PROJECTS + "/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed")
sys.path.insert(0, str(W / "code"))

PUBLISHED = 1500
CONST = "MIN_ARRAY_ELEMENT_SIZE"
SCAN = list(range(1000, 1751, 50))
WALK_LIMIT = 8                 # breakpoints to try per side before giving up
GUARD = (400, 3200)            # never test outside this, keeps runtime bounded

_recorded = []                 # (size, exon_fraction) for every cluster seen


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


def run_point(tag, value, record=False):
    base = W / "points" / tag
    if base.exists():
        shutil.rmtree(base)
    (base / "results").mkdir(parents=True)
    link = base / "data"
    if not link.exists():
        link.symlink_to(DATA)

    for m in [m for m in list(sys.modules) if m == "analyze_d4z4_arrays"]:
        del sys.modules[m]
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
    mod = importlib.import_module("analyze_d4z4_arrays")
    assert hasattr(mod, CONST), f"{CONST} is not a module attribute"
    setattr(mod, CONST, value)

    if record:
        _recorded.clear()
        original = mod.identify_arrays
        def recording_identify_arrays(clusters, gene_regions):
            for start, end in clusters or []:
                fr = mod.get_cluster_annotation_fractions(start, end, gene_regions)
                _recorded.append((end - start, fr["exon"]))
            return original(clusters, gene_regions)
        mod.identify_arrays = recording_identify_arrays

    os.environ["D4Z4_BASE_DIR"] = str(base)
    saved = sys.stdout
    try:
        with open(base / "stdout.log", "w") as fh:
            sys.stdout = fh
            mod.main()
    finally:
        sys.stdout = saved
    return read_results(base / "results")


def sig(agg, per_hap):
    return (tuple(sorted(agg.items())), tuple(sorted(per_hap.items())))


def diff_fields(a, b):
    return [k for k in a if a[k] != b[k]]


_cache = {}
def evaluate(value):
    if value not in _cache:
        agg, per_hap = run_point(f"{CONST}_{value}", value)
        _cache[value] = (agg, per_hap)
    return _cache[value]


def main():
    rows = []

    # -- Phase 1: baseline at the published value, recording cluster sizes ----
    print("=== phase 1: baseline at published value, recording cluster sizes ===", flush=True)
    base_agg, base_per_hap = run_point(f"{CONST}_{PUBLISHED}_rec", PUBLISHED, record=True)
    _cache[PUBLISHED] = (base_agg, base_per_hap)
    base_sig = sig(base_agg, base_per_hap)
    print(f"  baseline {base_agg}", flush=True)

    dense_thr = importlib.import_module("analyze_d4z4_arrays").EXON_DENSE_THRESHOLD
    qual = sorted({sz for sz, ex in _recorded if ex >= dense_thr})
    print(f"  clusters seen: {len(_recorded)}  qualifying (exon >= {dense_thr}): "
          f"{len(qual)} distinct sizes", flush=True)

    above = [s for s in qual if s >= PUBLISHED]
    below = [s for s in qual if s < PUBLISHED]
    s0 = above[0] if above else None
    t0 = below[-1] if below else None
    print(f"  smallest qualifying size >= {PUBLISHED}: {s0}", flush=True)
    print(f"  largest  qualifying size <  {PUBLISHED}: {t0}", flush=True)
    print(f"  => exon-dense SET is identical for T in "
          f"[{(t0 + 1) if t0 else GUARD[0]}, {s0 if s0 else GUARD[1]}]", flush=True)

    (W / "qualifying_cluster_sizes.json").write_text(json.dumps(
        {"exon_dense_threshold": dense_thr, "n_clusters_seen": len(_recorded),
         "n_qualifying_distinct_sizes": len(qual),
         "nearest_below_1500": below[-8:], "nearest_at_or_above_1500": above[:8]}, indent=2))

    # -- Phase 2: step-50 scan for the record --------------------------------
    print("\n=== phase 2: step-50 scan over [1000, 1750] ===", flush=True)
    for v in SCAN:
        agg, per_hap = evaluate(v)
        changed = sum(1 for h, n in per_hap.items() if base_per_hap.get(h) != n)
        same = sig(agg, per_hap) == base_sig
        rows.append(dict(phase="scan50", value=v, matches_published=int(same),
                         haps_changed=changed, **agg,
                         d_arrays=agg["arrays"] - base_agg["arrays"],
                         d_repeat_main=agg["repeat_main"] - base_agg["repeat_main"],
                         d_repeat_total=agg["repeat_total"] - base_agg["repeat_total"]))
        print(f"  {v:>5}  arrays={agg['arrays']:>5} repeat_main={agg['repeat_main']:>6} "
              f"repeat_total={agg['repeat_total']:>6} haps_changed={changed:>4} "
              f"{'SAME' if same else 'DIFFERS: ' + ','.join(diff_fields(base_agg, agg))}"
              + ("  <- published" if v == PUBLISHED else ""), flush=True)

    # -- Phase 3: walk the data breakpoints outward to the exact edges -------
    def walk(direction):
        """direction=+1 walks up, -1 walks down. Returns (edge, evidence rows)."""
        cands = above if direction > 0 else list(reversed(below))
        ev = []
        last_same = PUBLISHED
        for i, s in enumerate(cands[:WALK_LIMIT]):
            # T at which this breakpoint first alters the exon-dense set
            probe = s + 1 if direction > 0 else s
            if not (GUARD[0] <= probe <= GUARD[1]):
                print(f"  guard stop at probe={probe}", flush=True)
                break
            agg, per_hap = evaluate(probe)
            changed = sum(1 for h, n in per_hap.items() if base_per_hap.get(h) != n)
            same = sig(agg, per_hap) == base_sig
            fields = diff_fields(base_agg, agg)
            ev.append(dict(phase=f"walk{'up' if direction > 0 else 'down'}",
                           value=probe, matches_published=int(same),
                           haps_changed=changed, **agg,
                           d_arrays=agg["arrays"] - base_agg["arrays"],
                           d_repeat_main=agg["repeat_main"] - base_agg["repeat_main"],
                           d_repeat_total=agg["repeat_total"] - base_agg["repeat_total"]))
            print(f"  breakpoint size={s:<6} probe T={probe:<6} "
                  f"arrays={agg['arrays']:>5} repeat_main={agg['repeat_main']:>6} "
                  f"haps_changed={changed:>4} "
                  f"{'SAME' if same else 'DIFFERS: ' + ','.join(fields)}", flush=True)
            if not same:
                # plateau edge is one step inside this probe
                return (probe - 1 if direction > 0 else probe + 1), ev
            last_same = probe
        print(f"  no change within {WALK_LIMIT} breakpoints; edge is at least {last_same}",
              flush=True)
        return last_same, ev

    print("\n=== phase 3a: walk UP to the exact upper edge ===", flush=True)
    upper, ev_up = walk(+1)
    print("\n=== phase 3b: walk DOWN to the exact lower edge ===", flush=True)
    lower, ev_down = walk(-1)
    rows.extend(ev_up)
    rows.extend(ev_down)

    # -- Phase 4: confirm the edges are tight --------------------------------
    print("\n=== phase 4: confirm edges are tight ===", flush=True)
    for label, v in [("upper edge", upper), ("upper edge + 1", upper + 1),
                     ("lower edge", lower), ("lower edge - 1", lower - 1)]:
        if not (GUARD[0] <= v <= GUARD[1]):
            continue
        agg, per_hap = evaluate(v)
        same = sig(agg, per_hap) == base_sig
        print(f"  {label:<16} T={v:<6} {'SAME' if same else 'DIFFERS: ' + ','.join(diff_fields(base_agg, agg))}",
              flush=True)
        rows.append(dict(phase="confirm", value=v, matches_published=int(same),
                         haps_changed=sum(1 for h, n in per_hap.items() if base_per_hap.get(h) != n),
                         **agg,
                         d_arrays=agg["arrays"] - base_agg["arrays"],
                         d_repeat_main=agg["repeat_main"] - base_agg["repeat_main"],
                         d_repeat_total=agg["repeat_total"] - base_agg["repeat_total"]))

    print(f"\n=== RESULT ===", flush=True)
    print(f"  exact plateau: [{lower}, {upper}]   published = {PUBLISHED}", flush=True)
    print(f"  margin below published: {PUBLISHED - lower} bp", flush=True)
    print(f"  margin above published: {upper - PUBLISHED} bp", flush=True)
    print(f"  published value is {'INTERIOR to' if lower < PUBLISHED < upper else 'ON THE EDGE of'} the plateau",
          flush=True)

    cols = ["phase", "value", "matches_published", "n_hap", "arrays", "repeat_main",
            "repeat_total", "canonical", "degraded", "haps_changed",
            "d_arrays", "d_repeat_main", "d_repeat_total"]
    dest = W / "min_array_element_size_refined.tsv"
    dest.write_text("\t".join(cols) + "\n" +
                    "\n".join("\t".join(str(r[c]) for c in cols) for r in rows) + "\n")
    (W / "plateau.json").write_text(json.dumps(
        {"constant": CONST, "published": PUBLISHED, "plateau_lower": lower,
         "plateau_upper": upper, "margin_below": PUBLISHED - lower,
         "margin_above": upper - PUBLISHED,
         "published_is_interior": bool(lower < PUBLISHED < upper),
         "baseline": base_agg}, indent=2))
    print(f"\nwrote {dest} ({len(rows)} rows)", flush=True)


main()
