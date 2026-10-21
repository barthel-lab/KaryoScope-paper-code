#!/usr/bin/env python3
"""Print the A2 Step 3 per-sample values from the ROB configuration ($ROBS_CONFIG), one line per sample.

Usage: rob_samples.py FIELD [FIELD ...]
FIELD is a key of a2_step3.samples.<sample> (for example sst1_mid), or one of
  sample, scaffold, left, right, win_start, win_end
where scaffold, left (chrom1), right (chrom2) and the ISCN window come from the sample's breakpoint group.
Fields are separated by a space, samples in configuration order. The shell scripts here read their
per-sample values this way.
"""
import os
import sys

import yaml


def samples(path=None):
    """{sample: values} for a2_step3, with the breakpoint-group fields resolved."""
    c = yaml.safe_load(open(path or os.environ["ROBS_CONFIG"]))
    out = {}
    for s, d in c["a2_step3"]["samples"].items():
        g = c["breakpoint_groups"][d["breakpoint_group"]][0]
        out[s] = dict(d, sample=s, scaffold=g["scaffold"], left=g["chrom1"], right=g["chrom2"],
                      win_start=g["start"], win_end=g["end"])
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for s, d in samples().items():
        print(" ".join(str(d[f]) for f in sys.argv[1:]))
