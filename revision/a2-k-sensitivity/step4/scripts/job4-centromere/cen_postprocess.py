#!/usr/bin/env python3
"""Turn a k-sweep region annotation of the published centromere windows into the
representation Yi-An's pipeline consumed: 100 kb winner-take-all bins in the
contig frame (published: bin_features.py -b 100000, Feb 2026, plain mode),
trimmed per contig to [first centromere-class bin start, last centromere-class
bin start] exactly as pangenome_T2T_chromosome_centromere2.py did.

One adaptation to that script's get_simple_region(): v1's internal arm node was
'arm_multigroup1' (mapped to the q_arm class there); v2's is 'arm', which the
literal substring test would otherwise class as centromeric. 'arm' is mapped to
the arm class here. Labels are otherwise passed through as v2 names -- the
clustering is label-agnostic (features are whatever labels occur).
"""
import argparse, os, subprocess, sys, tempfile, collections

ap = argparse.ArgumentParser()
ap.add_argument("--bed", required=True, help="presmoothed region BED of the windows (window frame)")
ap.add_argument("--windows", required=True, help="common/cen_windows.tsv")
ap.add_argument("--chrom", required=True)
ap.add_argument("--bin-script", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--extents-out", required=True)
a = ap.parse_args()

def simple_region(name):
    n = name.lower()
    if "p_arm" in n: return "p_arm"
    if "q_arm" in n or "arm_multigroup1" in n or n == "arm": return "q_arm"
    if "tel" in n: return "telomere"
    return "centromere"

off = {}
for line in open(a.windows):
    p = line.rstrip("\n").split("\t")
    if p[1] == a.chrom and p[8] == "1":
        off[p[0]] = int(p[2])

tmp = tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed", dir=os.path.dirname(a.out))
with open(a.bed) as fh:
    for line in fh:
        c, s, e, lab = line.rstrip("\n").split("\t")[:4]
        o = off[c]
        tmp.write(f"{c}\t{int(s)+o}\t{int(e)+o}\t{lab}\n")
tmp.close()
binned = tmp.name + ".b100k"
subprocess.run([sys.executable, a.bin_script, "-b", "100000", "-i", tmp.name, "-o", binned], check=True)

rows = collections.defaultdict(list)
with open(binned) as fh:
    for line in fh:
        c, s, e, lab = line.rstrip("\n").split("\t")[:4]
        rows[c].append((int(s), int(e), lab))
n_out = 0
with open(a.out, "w") as o, open(a.extents_out, "w") as ex:
    ex.write("contig\tchrom\twindow_start\twindow_end\ttrim_start\ttrim_end\tn_bins_kept\n")
    for c in off:
        r = rows.get(c, [])
        cen = [s for s, e, lab in r if simple_region(lab) == "centromere"]
        if not cen:
            ex.write(f"{c}\t{a.chrom}\t{off[c]}\t{r[-1][1] if r else off[c]}\tNA\tNA\t0\n"); continue
        lo, hi = min(cen), max(cen)
        kept = [(s, e, lab) for s, e, lab in r if lo <= s <= hi]
        for s, e, lab in kept:
            o.write(f"{c}\t{s}\t{e}\t{lab}\n"); n_out += 1
        ex.write(f"{c}\t{a.chrom}\t{off[c]}\t{r[-1][1]}\t{lo}\t{kept[-1][1]}\t{len(kept)}\n")
os.unlink(tmp.name); os.unlink(binned)
print(f"{a.chrom}: {len(off)} windows, {n_out} binned rows -> {a.out}", flush=True)
