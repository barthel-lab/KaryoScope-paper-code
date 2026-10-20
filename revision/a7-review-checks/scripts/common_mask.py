"""Shared helpers: hierarchy parsing and the A4 N-mask convention.

N-mask convention (copied from a4_smoothing/run_a4_novel_v3.sh and run_a4_composition.sh):
every N run in the query FASTA is extended 30 bp LEFT (fill_gaps' k-1 at k=31) and the
overlap of each BED interval with those masked intervals is subtracted. Sequences the
database excludes (chrM) are dropped entirely.
"""
import collections, gzip, re

def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def parse_hierarchy(path):
    """Return {feature_set: (nodes, leaves, parent_of)}."""
    nodes = collections.defaultdict(set); parents = collections.defaultdict(set)
    parent_of = collections.defaultdict(dict)
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3 or p[0] == "feature_set": continue
        fs, child, parent = p[0], p[1], p[2]
        nodes[fs].add(child); nodes[fs].add(parent); parents[fs].add(parent)
        parent_of[fs][child] = parent
    return {fs: (nodes[fs], nodes[fs] - parents[fs], parent_of[fs]) for fs in nodes}

def make_level(nodes, leaves):
    def lvl(lab):
        if lab.lower() == "novel": return "novel"
        if lab in leaves: return "leaf"
        if lab in nodes: return "internal"
        return "other"
    return lvl

def ancestors(parent_of, node):
    out = [node]
    while node in parent_of:
        node = parent_of[node]; out.append(node)
    return out

def lca(parent_of, a, b):
    anc_b = set(ancestors(parent_of, b))
    for x in ancestors(parent_of, a):
        if x in anc_b: return x
    return None

def n_runs_from_fasta(path, ext_left=30):
    """{seq: [(start-ext_left, end), ...]} for every N run (case-insensitive)."""
    mask = collections.defaultdict(list); name = None; pos = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]; pos = 0; continue
            s = line.strip()
            for m in re.finditer(r"N+", s, flags=re.I):
                mask[name].append((max(0, pos + m.start() - ext_left), pos + m.end()))
            pos += len(s)
    # merge touching/overlapping
    out = {}
    for c, iv in mask.items():
        iv.sort(); merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else: merged.append((a, b))
        out[c] = merged
    return out

def n_runs_from_bed(path, ext_left=30):
    mask = collections.defaultdict(list)
    for line in open(path):
        c, s, e = line.split("\t")[:3]
        mask[c].append((max(0, int(s) - ext_left), int(e)))
    return dict(mask)

def overlap(mask_list, s, e):
    tot = 0
    for a, b in mask_list:
        lo, hi = max(s, a), min(e, b)
        if hi > lo: tot += hi - lo
    return tot
