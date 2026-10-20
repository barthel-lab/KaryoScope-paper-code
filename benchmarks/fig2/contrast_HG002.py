#!/usr/bin/env python3
"""
contrast_HG002.py

Compare truth annotations (RepeatMasker or CenSat) against KaryoScope
predictions and classify each region's compatibility using the feature
hierarchy.

Compatibility categories:
  - equivalent:   truth and prediction are the same label
  - compatible:   truth is a descendant of (more specific than) prediction
  - novel:        KaryoScope predicted "novel"

For incompatible predictions, the category is derived from the
top-level "bucket" (first child of 'categorized') each label falls
under in the hierarchy tree:

  Same bucket:
    wrong_{bucket}           e.g., wrong_centromeric, wrong_repeat

  Different bucket:
    {truth_bucket}_as_{pred_bucket}
                             e.g., centromeric_as_arm, repeat_as_nonrepeat

This generalizes across feature sets (region, repeat, etc.) without
hardcoding any particular tree structure.

Label merging (--merge):
  Multiple labels can be merged into a single label for evaluation
  purposes.  For example, Simple_repeat, Satellite, Low_complexity,
  and Noninterspersed can be merged into "tandem_repeat" when these
  categories share too much k-mer content to be reliably distinguished.

  Both truth and prediction labels are remapped before compatibility
  is assessed.  The merged label inherits the intersection of its
  constituents' ancestors (the most conservative correct placement
  in the hierarchy).

Usage:
    python contrast_HG002.py --input comparison.bed \\
        --hierarchy hierarchy.txt --feature-set region --output contrast.tsv

    python contrast_HG002.py --input comparison.bed \\
        --hierarchy hierarchy.txt --feature-set repeat \\
        --merge "Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat" \\
        --output contrast.tsv
"""

import sys
import argparse
from collections import defaultdict


def parse_hierarchy(hierarchy_file, feature_set):
    """
    Parse the KaryoScope hierarchy file for a given feature set.

    Returns:
      child_to_parent: {child: parent} mapping
      ancestors_map:   {node: set of all ancestors including self}
      bucket_map:      {node: top-level child of 'categorized' it falls under}
    """
    child_to_parent = {}
    with open(hierarchy_file, "r") as f:
        for line in f:
            # Skip header
            if line.startswith("feature_set"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            fs, child, parent = parts[0], parts[1], parts[2]
            if fs == feature_set:
                child_to_parent[child] = parent

    # Build transitive closure: for each node, find all ancestors
    ancestors_cache = {}

    def get_ancestors(node):
        if node in ancestors_cache:
            return ancestors_cache[node]
        ancestors = {node}
        if node in child_to_parent:
            parent = child_to_parent[node]
            ancestors.update(get_ancestors(parent))
        ancestors_cache[node] = ancestors
        return ancestors

    all_nodes = set(child_to_parent.keys()) | set(child_to_parent.values())
    ancestors_map = {}
    for node in all_nodes:
        ancestors_map[node] = get_ancestors(node)

    # Build bucket map: for each node, find its top-level ancestor
    # (the child of 'categorized' it falls under)
    bucket_map = {}
    for node in all_nodes:
        if node == "categorized":
            continue
        # Walk up until we find a node whose parent is 'categorized'
        current = node
        while current in child_to_parent and child_to_parent[current] != "categorized":
            current = child_to_parent[current]
        bucket_map[node] = current

    return child_to_parent, ancestors_map, bucket_map


def parse_merge_specs(merge_args):
    """
    Parse --merge arguments into a remap dictionary.

    Each argument has the format "label1,label2,...:merged_name".

    Returns:
      remap:        {original_label: merged_name} for all constituent labels
      merge_groups: {merged_name: [list of original labels]}
    """
    remap = {}
    merge_groups = {}

    if not merge_args:
        return remap, merge_groups

    for spec in merge_args:
        if ":" not in spec:
            print(f"Error: Invalid --merge format '{spec}'. "
                  f"Expected 'label1,label2:merged_name'.", file=sys.stderr)
            sys.exit(1)

        labels_str, merged_name = spec.rsplit(":", 1)
        labels = [l.strip() for l in labels_str.split(",")]

        if len(labels) < 2:
            print(f"Error: --merge requires at least 2 labels, "
                  f"got '{spec}'.", file=sys.stderr)
            sys.exit(1)

        for label in labels:
            if label in remap:
                print(f"Error: Label '{label}' appears in multiple "
                      f"--merge groups.", file=sys.stderr)
                sys.exit(1)
            remap[label] = merged_name

        merge_groups[merged_name] = labels

    return remap, merge_groups


def apply_merges(ancestors_map, bucket_map, merge_groups):
    """
    Add merged labels to ancestors_map and bucket_map.

    The merged label's ancestors are the intersection of its
    constituents' ancestor sets (excluding the constituents themselves),
    plus the merged label itself.  This is the most conservative
    placement: the merged label is compatible with any ancestor shared
    by all constituents.

    The merged label's bucket is the shared bucket of all constituents.
    If they don't share a bucket, an error is raised.
    """
    for merged_name, labels in merge_groups.items():
        # Compute ancestor intersection (excluding constituent labels)
        constituent_set = set(labels)
        ancestor_sets = []
        for label in labels:
            if label in ancestors_map:
                ancestor_sets.append(ancestors_map[label] - constituent_set)
            else:
                print(f"Warning: Merge label '{label}' not found in "
                      f"hierarchy; skipping its ancestors.",
                      file=sys.stderr)

        if ancestor_sets:
            shared_ancestors = set.intersection(*ancestor_sets)
        else:
            shared_ancestors = set()

        # Add merged label to ancestors map
        ancestors_map[merged_name] = shared_ancestors | {merged_name}

        # Verify shared bucket
        buckets = set()
        for label in labels:
            if label in bucket_map:
                buckets.add(bucket_map[label])

        if len(buckets) > 1:
            print(f"Error: --merge labels for '{merged_name}' span "
                  f"multiple buckets: {buckets}. Cannot merge across "
                  f"top-level categories.", file=sys.stderr)
            sys.exit(1)

        if buckets:
            bucket_map[merged_name] = buckets.pop()

    return ancestors_map, bucket_map


def get_compatibility(truth, prediction, ancestors_map, bucket_map):
    """
    Classify the compatibility between a truth annotation and a
    KaryoScope prediction.
    """
    if prediction == "novel":
        return "novel"

    if truth == prediction:
        return "equivalent"

    # Check ancestor-based compatibility: truth is more specific than prediction
    if truth in ancestors_map and prediction in ancestors_map.get(truth, set()):
        return "compatible"

    # Incompatible: classify by bucket
    truth_bucket = bucket_map.get(truth, truth)
    pred_bucket = bucket_map.get(prediction, prediction)

    if truth_bucket == pred_bucket:
        return f"wrong_{truth_bucket}"
    else:
        return f"{truth_bucket}_as_{pred_bucket}"


def main():
    parser = argparse.ArgumentParser(
        description="Compare truth annotations against KaryoScope predictions "
                    "using the feature hierarchy."
    )
    parser.add_argument(
        "--input", required=True,
        help="Input comparison BED file (6 columns: chrom, start, stop, "
             "truth, presmoothed, smoothed)."
    )
    parser.add_argument(
        "--hierarchy", required=True,
        help="KaryoScope hierarchy file (feature_set, child, parent)."
    )
    parser.add_argument(
        "--feature-set", required=True,
        help="Feature set to use from the hierarchy file "
             "(e.g., 'region', 'repeat')."
    )
    parser.add_argument(
        "--merge", action="append", default=None,
        help="Merge labels for evaluation: 'label1,label2,...:merged_name'. "
             "Can be specified multiple times. Both truth and prediction "
             "labels are remapped before compatibility assessment."
    )
    parser.add_argument(
        "--output", required=True,
        help="Output contrast summary file."
    )
    args = parser.parse_args()

    _, ancestors_map, bucket_map = parse_hierarchy(
        args.hierarchy, args.feature_set
    )

    # Parse and apply label merges
    remap, merge_groups = parse_merge_specs(args.merge)
    if merge_groups:
        ancestors_map, bucket_map = apply_merges(
            ancestors_map, bucket_map, merge_groups
        )
        merged_labels = ", ".join(
            f"{k} <- [{', '.join(v)}]" for k, v in merge_groups.items()
        )
        print(f"Label merges applied: {merged_labels}", file=sys.stderr)

    pre_features_to_length = defaultdict(int)
    smoothed_features_to_length = defaultdict(int)

    with open(args.input, "r") as input_file:
        for line in input_file:
            parts = line.strip().split()
            if len(parts) < 6:
                print(f"Warning: Skipping malformed line: {line.strip()}",
                      file=sys.stderr)
                continue
            chrom, start, stop = parts[0], int(parts[1]), int(parts[2])
            truth, presmoothed, smoothed = parts[3], parts[4], parts[5]
            length = stop - start

            # Skip rows where any column is '0' — this is a bedtools
            # unionbedg artifact from bases at contig ends that are not
            # covered by one of the input BED files (e.g., KaryoScope
            # cannot annotate the last k-1 bases of each contig).
            if "0" in (truth, presmoothed, smoothed):
                continue

            # Skip regions marked as 'exclude' in the truth annotations
            # (e.g., chrM, which KaryoScope does not meaningfully annotate).
            if truth == "exclude":
                continue

            # Apply label remapping for merged categories
            truth = remap.get(truth, truth)
            presmoothed = remap.get(presmoothed, presmoothed)
            smoothed = remap.get(smoothed, smoothed)

            pre_compat = get_compatibility(
                truth, presmoothed, ancestors_map, bucket_map
            )
            smoothed_compat = get_compatibility(
                truth, smoothed, ancestors_map, bucket_map
            )

            pre_features_to_length[
                (truth, presmoothed, pre_compat)
            ] += length
            smoothed_features_to_length[
                (truth, smoothed, smoothed_compat)
            ] += length

    with open(args.output, "w") as output_file:
        # Header
        output_file.write(
            "stage\ttruth\tprediction\tlength\tcompatibility\n"
        )
        for (truth, prediction, compat), length in sorted(
            pre_features_to_length.items(),
            key=lambda x: x[1], reverse=True
        ):
            output_file.write(
                f"presmoothed\t{truth}\t{prediction}\t{length}\t{compat}\n"
            )
        for (truth, prediction, compat), length in sorted(
            smoothed_features_to_length.items(),
            key=lambda x: x[1], reverse=True
        ):
            output_file.write(
                f"smoothed\t{truth}\t{prediction}\t{length}\t{compat}\n"
            )

    print("Contrast complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
