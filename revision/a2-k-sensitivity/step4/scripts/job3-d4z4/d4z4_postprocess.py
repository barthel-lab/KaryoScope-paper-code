#!/usr/bin/env python3
"""Turn a k-sweep annotation of the 500 kb q-terminal windows into the exact
input representation the published D4Z4 caller consumed.

Published chain (KaryoScope-archive, Feb 2026): whole-assembly presmoothed BED,
scaffolded -> bin_features.py -b 100 (plain max-overlap, contig frame) ->
bedtools intersect with the [L-500000, L) mask -> normalize_and_concat.py
(subtract L-500000, prefix chr<N>_hap<H>_).  Reproduced here on the window
annotation by shifting into the contig frame BEFORE binning, so bin boundaries
land where the published ones did, then clipping and normalising back.

The six-k family carries v2 label names; the caller string-matches v1 names.
Only the labels the caller tests are renamed (a 1:1 rename, nothing merged):
  gene:         exon/intron/intergenic -> *_specific
  region:       bSat -> bsat_specific, ct -> ct_specific
  subtelomeric: TAR1 -> TAR1_specific
  repeat:       LINE/LTR/SINE -> *_specific
Everything else, including 'novel', is passed through unchanged.
"""
import argparse, os, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--annot-dir", required=True)
ap.add_argument("--db", required=True, help="database id in the BED filenames")
ap.add_argument("--region-annot-dir", default=None, help="optional: take the REGION set from this dir instead")
ap.add_argument("--region-db", default=None, help="db id in the region BED filenames (with --region-annot-dir)")
ap.add_argument("--windows", required=True, help="common/d4z4_windows.tsv")
ap.add_argument("--bin-script", required=True)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()

RENAME = {
    "gene": {"exon": "exon_specific", "intron": "intron_specific", "intergenic": "intergenic_specific"},
    "region": {"bSat": "bsat_specific", "ct": "ct_specific", "bsat": "bsat_specific", "p_arm": "p_arm_specific", "q_arm": "q_arm_specific",
               "censat": "censat_specific", "gsat": "gsat_specific", "hsat1A": "hsat1A_specific", "hsat1B": "hsat1B_specific",
               "hsat2": "hsat2_specific", "hsat3": "hsat3_specific", "active": "active_specific", "inactive": "inactive_specific",
               "divergent": "divergent_specific", "monomeric": "monomeric_specific", "noncentromeric": "noncentromeric_specific",
               "arm": "arm_multigroup1", "hor": "hor_multigroup1", "asat": "asat_multigroup1", "hsat": "hsat_multigroup1",
               "hsat1": "hsat1_multigroup1", "centromeric": "centromeric_multigroup1"},
    "subtelomeric": {"TAR1": "TAR1_specific"},
    "repeat": {"LINE": "LINE_specific", "LTR": "LTR_specific", "SINE": "SINE_specific"},
}
WIN = 500000
# published_name, pansn, chrom, hap, fasta, L, start, end
info = {}
for line in open(a.windows):
    p = line.rstrip("\n").split("\t")
    info[p[1]] = dict(pub=p[0], chrom=p[2], L=int(p[5]), off=int(p[6]))

os.makedirs(a.outdir, exist_ok=True)
for chrom in ("chr4", "chr10"):
    for fs in ("region", "subtelomeric", "gene", "repeat"):
        if fs == "region" and a.region_annot_dir:
            src = os.path.join(a.region_annot_dir, f"windows_{chrom}.{a.region_db}.{fs}.presmoothed.bed")
        else:
            src = os.path.join(a.annot_dir, f"windows_{chrom}.{a.db}.{fs}.presmoothed.bed")
        if not os.path.exists(src):
            sys.exit(f"missing {src}")
        tmp = tempfile.NamedTemporaryFile("w", delete=False, suffix=".bed", dir=a.outdir)
        with open(src) as fh:
            for line in fh:
                c, s, e, lab = line.rstrip("\n").split("\t")[:4]
                off = info[c]["off"]
                tmp.write(f"{c}\t{int(s)+off}\t{int(e)+off}\t{lab}\n")
        tmp.close()
        binned = tmp.name + ".b100"
        subprocess.run([sys.executable, a.bin_script, "-b", "100", "-i", tmp.name, "-o", binned], check=True)
        dest = os.path.join(a.outdir, f"pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.{fs}.bed")
        n = 0
        with open(binned) as fh, open(dest, "w") as o:
            for line in fh:
                c, s, e, lab = line.rstrip("\n").split("\t")[:4]
                off = info[c]["off"]
                s = int(s) - off; e = int(e) - off
                if e <= 0: continue
                s = max(0, s); e = min(WIN, e)   # bedtools-intersect clipping to the mask
                if e <= s: continue
                lab = RENAME[fs].get(lab, lab)
                o.write(f"{info[c]['pub']}\t{s}\t{e}\t{lab}\n"); n += 1
        os.unlink(tmp.name); os.unlink(binned)
        print(f"{chrom} {fs}: {n} rows -> {dest}", flush=True)
