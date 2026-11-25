"""Display names of KaryoScope feature labels: what a figure prints for a label.

The labels in the databases and in every data file stay as they are; a figure script that is given
--display-names passes each label it draws through display_name() instead of printing it raw. One
table serves every figure:

    from display_names import display_name, MULTIGROUP
    display_name("active_hor")                  # 'active HOR'   (v2 database)
    display_name("active_specific", "v1")       # 'active HOR'   (v1 database)
    display_name("asat_multigroup1", "v1")      # 'aSat'
    display_name("telomere_like")               # 'telomere-like'

`vocabulary` is "v2" (KS_human_CHM13_v2: bare names, e.g. active_hor, hor, mon) or "v1"
(KS_human_CHM13: <name>_specific leaves and <name>_multigroup1 parents, e.g. active_specific,
inactive_specific, hor_multigroup1). The two differ only in the region feature set, where "hor" is
the inactive-HOR leaf in v2 but the alpha-satellite HOR parent in v1.

The v1 and v2 tables agree with karyoscope_feature_name_mapping.tsv (here; the v1 -> v2 label table): every
pair in it gets one display name. v1 noncentromeric_specific is v2 rDNA (the rDNA intervals the v1 region
source dropped) and prints "rDNA".

Legend order (display_sort_key): the region feature set in the order of REGION_ORDER, which is the order of
the Fig. 1 colour key -- the leaves (active HOR, divergent HOR, inactive HOR, mixed aSat, monomeric, bSat,
cenSat, ct, gSat, HSat1A, HSat1B, HSat2, HSat3, rDNA, p arm, q arm), then the parents (centromeric, aSat,
aSat HOR, HSat, HSat1, arm); any other name after these, alphabetically.

display_text(phrase) rewrites a printed phrase: each token that is a feature label (with or without its
_specific / _multigroup1 suffix) becomes its display name, and the Greek forms beta-satellite / beta-sat become
bSat; every other word is kept.

shared_parent_name(labels) names a legend swatch that several parents share: the deepest of them in the v1
hierarchy (vendored/KaryoScope-archive-resources/resources/CHM13/CHM13_*.priority.txt), the root excluded --
the one closest to the leaves drawn (acrocentric, not autosome, for chr13/chr14/chr21).

A parent label (<name>_multigroup1) prints as the parent's display name alone: a parent name already
means "not resolved to a more specific label". MULTIGROUP is the one text a figure may still print
after a parent row's k-mer count (the Fig. 1b tree, by option). A label not in the tables prints
with its underscores as spaces, and is recorded in UNLISTED so a figure script can report it.
"""

MULTIGROUP = "(multigroup)"

# region feature set, v2 database (bare names)
REGION_V2 = {
    "active_hor": "active HOR",
    "hor": "inactive HOR",
    "dhor": "divergent HOR",
    "mon": "monomeric",
    "mixedAlpha": "mixed aSat",
    "aSat": "aSat",
    "alpha_hor": "aSat HOR",
    "bSat": "bSat",
    "gSat": "gSat",
    "cenSat": "cenSat",
    "ct": "ct",
    "HSat": "HSat", "HSat1": "HSat1", "HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2", "HSat3": "HSat3",
    "rDNA": "rDNA",
    "p_arm": "p arm", "q_arm": "q arm", "arm": "arm",
    "centromeric": "centromeric",
}

# region feature set, v1 database (the names without _specific / _multigroup1)
REGION_V1 = {
    "active": "active HOR",
    "inactive": "inactive HOR",
    "divergent": "divergent HOR",
    "monomeric": "monomeric",
    "asat": "aSat",
    "hor": "aSat HOR",                 # a parent in v1 (hor_multigroup1)
    "bsat": "bSat",
    "gsat": "gSat",
    "censat": "cenSat",
    "ct": "ct",
    "hsat": "HSat", "hsat1": "HSat1", "hsat1A": "HSat1A", "hsat1B": "HSat1B", "hsat2": "HSat2", "hsat3": "HSat3",
    "p_arm": "p arm", "q_arm": "q arm", "arm": "arm",
    "centromeric": "centromeric",
    "noncentromeric": "rDNA",              # v2 rDNA (karyoscope_feature_name_mapping.tsv)
}

# every other feature set (the same in v1 and v2)
OTHER = {
    "canonical_telomere": "canonical telomere",
    "noncanonical_telomere": "noncanonical telomere",
    "telomere_like": "telomere-like",
    "nonsubtelomeric": "nonsubtelomeric",
    "nonrepeat": "nonrepeat",
    "nonacrocentric": "nonacrocentric",
    "Class_I_Retrotransposition": "Class I Retrotransposition",
    "Class_II_DNA_Transposition": "Class II DNA Transposition",
    "Transposable_Element": "Transposable Element",
    "Interspersed_Repeat": "Interspersed Repeat",
    "LINE-dependent_Retroposon": "LINE-dependent Retroposon",
    "tandem_repeat": "tandem repeat",
    "Simple_repeat": "Simple repeat",
    "Low_complexity": "Low complexity",
    **{k: k for k in ("LINE", "SINE", "LTR", "DNA", "RC", "Retroposon", "rRNA", "snRNA", "scRNA", "srpRNA",
                      "tRNA", "Unknown", "exon", "intron", "intergenic", "DJ", "PJ", "PHR", "SST1", "ITS",
                      "TAR1", "novel")},
}

UNLISTED = set()


def display_name(label, vocabulary="v2"):
    """The name a figure prints for a feature label (see the module docstring)."""
    if vocabulary not in ("v1", "v2"):
        raise ValueError(f"vocabulary must be 'v1' or 'v2', not {vocabulary!r}")
    base = label
    if base.endswith("_multigroup1"):
        base = base[: -len("_multigroup1")]
    elif base.endswith("_specific"):
        base = base[: -len("_specific")]
    region = REGION_V1 if vocabulary == "v1" else REGION_V2
    if base in region:
        name = region[base]
    elif base in OTHER:
        name = OTHER[base]
    else:
        UNLISTED.add(label)
        name = base.replace("_", " ")
    return name


REGION_ORDER = ["active HOR", "divergent HOR", "inactive HOR", "mixed aSat", "monomeric", "bSat", "cenSat", "ct",
                "gSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "rDNA", "p arm", "q arm",
                "centromeric", "aSat", "aSat HOR", "HSat", "HSat1", "arm"]


def display_sort_key(name):
    """Sort key of a printed (display) name: REGION_ORDER first, then the rest alphabetically."""
    return (REGION_ORDER.index(name), "") if name in REGION_ORDER else (len(REGION_ORDER), name.lower())


GREEK = {"β-satellite": "bSat", "β-sat": "bSat"}


def display_text(phrase, vocabulary="v1"):
    """A printed phrase with its feature labels as display names (see the module docstring)."""
    import re
    for g, d in GREEK.items():
        phrase = phrase.replace(g, d)
    region = REGION_V1 if vocabulary == "v1" else REGION_V2

    def sub(m):
        tok = m.group(0)
        base = re.sub(r"_(specific|multigroup1)$", "", tok)
        if base in region or base in OTHER:
            return display_name(tok, vocabulary)
        return tok
    return re.sub(r"[A-Za-z][A-Za-z0-9_]*", sub, phrase)


_PARENTS = None


def _v1_parents():
    global _PARENTS
    if _PARENTS is None:
        from pathlib import Path
        _PARENTS = {}
        d = Path(__file__).resolve().parent.parent / "vendored/KaryoScope-archive-resources/resources/CHM13"
        for f in sorted(d.glob("CHM13_*.priority.txt")):
            for line in f.read_text().splitlines():
                p = line.split("\t")
                if len(p) >= 3:
                    _PARENTS[p[0]] = p[2]
    return _PARENTS


def shared_parent_name(labels, vocabulary="v1"):
    """The display name for a swatch several parent labels share: the deepest of them, root excluded."""
    import re
    parents = _v1_parents()

    def depth(n):
        k = 0
        while n in parents and k < 50:
            n, k = parents[n], k + 1
        return k if n == "categorized" else -1
    names = sorted({re.sub(r"_multigroup\d*$", "", l) for l in labels if l != "categorized"})
    if not names:
        return None
    best = max(names, key=lambda n: (depth(n), [-ord(c) for c in n]))
    return display_name(best, vocabulary)
