#!/usr/bin/env python3
"""Inject controlled SNP/indel divergence into a FASTA, with truth projection (A5).

Writes, for --out PREFIX:
  PREFIX.fa                       mutated FASTA
  PREFIX.summary.tsv              per-sequence counts: subs, ins, del, bases inserted/deleted
  PREFIX.inserted.bed             inserted-base intervals in MUTATED coordinates
                                  (exclude these bases from concordance scoring; they
                                  have no original label)
  PREFIX.<name>.projected.bed     each --bed projected into MUTATED coordinates, so the
                                  original annotation can be scored per base against the
                                  mutated copy's annotation

Deterministic for a given --seed. Substitutions and indels are drawn independently per
1 Mb chunk (binomial counts, uniform positions), so rates are per-base probabilities.
Indel lengths are geometric (mean ~1/p) capped at --indel-max; insertions and deletions
are equally likely. Non-ACGT bases (N runs) are never substituted; indels are not placed
in them either, so gap structure is preserved.

Example (substitution-only series):
  python3 mutate_fasta.py -i chm13.fa.gz -o chm13_sub0.5 --sub-rate 0.005 --indel-rate 0 \
      --seed 42 --bed repeat=chm13_repeat.bed --bed region=chm13_region.bed
"""
import argparse
import bisect
import gzip
import sys

import numpy as np

CHUNK = 1_000_000
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
# For each base byte, the 3 alternatives to substitute with (uppercase handled; the
# sequence is uppercased on read).
ALT = {b: [x for x in b"ACGT" if x != b] for b in b"ACGT"}


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, parts = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def geometric_lengths(rng, n, p, cap):
    if n == 0:
        return np.empty(0, dtype=np.int64)
    return np.minimum(rng.geometric(p, size=n), cap)


def mutate_sequence(seq, rng, sub_rate, indel_rate, indel_p, indel_max, out_fh,
                    inserted_bed, chrom):
    """Mutate one sequence, streaming the result to out_fh (already '>name'-headered).

    Returns (breakpoints, stats). breakpoints is a list of
    (orig_start, orig_end, mut_start, cum_delta_after) per indel, in original-coordinate
    order: a deletion consumes [orig_start, orig_end) and emits nothing; an insertion has
    orig_start == orig_end and emits bases at mut_start.
    """
    raw = np.frombuffer(seq.encode(), dtype=np.uint8)
    n = raw.size
    breakpoints = []
    stats = {"sub": 0, "ins": 0, "del": 0, "ins_bp": 0, "del_bp": 0}
    mut_pos = 0          # bases written so far (mutated coordinate)
    cum_delta = 0
    line_buf = bytearray()

    def emit(arr):
        nonlocal mut_pos
        line_buf.extend(arr.tobytes() if isinstance(arr, np.ndarray) else arr)
        mut_pos += len(arr)
        while len(line_buf) >= 60:
            out_fh.write(line_buf[:60].decode())
            out_fh.write("\n")
            del line_buf[:60]

    for cstart in range(0, n, CHUNK):
        chunk = raw[cstart:cstart + CHUNK].copy()
        clen = chunk.size
        acgt = np.isin(chunk, BASES)
        n_eligible = int(acgt.sum())
        if n_eligible == 0:
            emit(chunk)
            continue
        eligible_idx = np.flatnonzero(acgt)

        # --- substitutions ---
        n_sub = rng.binomial(n_eligible, sub_rate) if sub_rate > 0 else 0
        if n_sub:
            pos = rng.choice(eligible_idx, size=n_sub, replace=False)
            for p in pos:
                chunk[p] = ALT[chunk[p]][rng.integers(3)]
            stats["sub"] += int(n_sub)

        # --- indels ---
        n_ind = rng.binomial(n_eligible, indel_rate) if indel_rate > 0 else 0
        if n_ind == 0:
            emit(chunk)
            continue
        ind_pos = np.sort(rng.choice(eligible_idx, size=n_ind, replace=False))
        is_ins = rng.random(n_ind) < 0.5
        lens = geometric_lengths(rng, n_ind, indel_p, indel_max)

        cursor = 0
        for p, ins, L in zip(ind_pos.tolist(), is_ins.tolist(), lens.tolist()):
            if p < cursor:            # swallowed by a previous deletion
                continue
            emit(chunk[cursor:p])
            L = int(L)
            if ins:
                orig = cstart + p
                mstart = mut_pos
                emit(BASES[rng.integers(0, 4, size=L)])
                cum_delta += L
                breakpoints.append((orig, orig, mstart, cum_delta))
                inserted_bed.write(f"{chrom}\t{mstart}\t{mstart + L}\n")
                stats["ins"] += 1
                stats["ins_bp"] += L
                cursor = p
            else:
                L = min(L, clen - p)  # cap at chunk end
                orig = cstart + p
                cum_delta -= L
                breakpoints.append((orig, orig + L, mut_pos, cum_delta))
                stats["del"] += 1
                stats["del_bp"] += L
                cursor = p + L
        emit(chunk[cursor:])

    if line_buf:
        out_fh.write(line_buf.decode())
        out_fh.write("\n")
    return breakpoints, stats


def project_point(bps_starts, bps, p):
    """Original coordinate -> mutated coordinate, clamping inside deletions."""
    i = bisect.bisect_right(bps_starts, p) - 1
    if i < 0:
        return p
    orig_start, orig_end, mut_start, cum_after = bps[i]
    if p < orig_end:                  # inside a deleted segment -> clamp to its start
        return mut_start
    return p + cum_after


def project_bed(in_path, out_path, bp_map):
    op = gzip.open if in_path.endswith(".gz") else open
    with op(in_path, "rt") as fh, open(out_path, "w") as out:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            bps = bp_map.get(chrom)
            if bps is None:
                out.write(line)       # chromosome untouched (or absent): pass through
                continue
            starts, tuples = bps
            ms, me = project_point(starts, tuples, s), project_point(starts, tuples, e)
            if me > ms:
                out.write("\t".join([chrom, str(ms), str(me)] + f[3:]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="FASTA (.gz ok)")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--sub-rate", type=float, required=True,
                    help="per-base substitution probability, e.g. 0.005 for 0.5%%")
    ap.add_argument("--indel-rate", type=float, required=True,
                    help="per-base indel-event probability (0 for substitution-only)")
    ap.add_argument("--indel-p", type=float, default=1 / 3,
                    help="geometric length parameter (mean length ~1/p; default 1/3)")
    ap.add_argument("--indel-max", type=int, default=50, help="max indel length")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--bed", action="append", default=[], metavar="NAME=PATH",
                    help="truth BED to project into mutated coordinates (repeatable)")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    bp_map = {}
    with open(f"{args.out}.fa", "w") as fa, \
         open(f"{args.out}.inserted.bed", "w") as ins_bed, \
         open(f"{args.out}.summary.tsv", "w") as summ:
        summ.write("chrom\tlength\tsubs\tins\tdel\tins_bp\tdel_bp\tmut_length\n")
        for name, seq in read_fasta(args.input):
            fa.write(f">{name}\n")
            bps, st = mutate_sequence(seq, rng, args.sub_rate, args.indel_rate,
                                      args.indel_p, args.indel_max, fa, ins_bed, name)
            bp_map[name] = ([b[0] for b in bps], bps)
            mut_len = len(seq) + st["ins_bp"] - st["del_bp"]
            summ.write(f"{name}\t{len(seq)}\t{st['sub']}\t{st['ins']}\t{st['del']}\t"
                       f"{st['ins_bp']}\t{st['del_bp']}\t{mut_len}\n")
            print(f"{name}: {st}", file=sys.stderr)

    for spec in args.bed:
        bed_name, path = spec.split("=", 1)
        project_bed(path, f"{args.out}.{bed_name}.projected.bed", bp_map)


if __name__ == "__main__":
    main()
