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
pair in it gets one display name, except v1 noncentromeric_specific, which the table maps to v2 rDNA (the
rDNA intervals the v1 region source dropped) and which prints "noncentromeric" here.

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
    "noncentromeric": "noncentromeric",
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
