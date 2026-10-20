# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/tsv_generator.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Generate TSV files from YAML breakpoint/label definitions.

The KaryoScope scripts (KaryoScope_assembly_contig_zoom_plot.py,
KaryoScope_rearrangement_plot.py) expect regions and labels as TSV files.
This module converts the YAML-defined breakpoints and labels into the
exact TSV formats those scripts consume.
"""

import os


REGIONS_HEADER = (
    "scaffold\tregion_type\tcenter\tstart\tend\t"
    "chrom1\tchrom2\tarm1\tarm2\tiscn_notation\toriginal_scaffold"
)


def write_regions_tsv(breakpoints: list, output_path: str) -> str:
    """Write a regions TSV from a list of breakpoint dicts.

    Output format (11 columns, tab-separated, with header):
        scaffold  region_type  center  start  end  chrom1  chrom2
        arm1  arm2  iscn_notation  original_scaffold
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(REGIONS_HEADER + "\n")
        for bp in breakpoints:
            f.write(
                f"{bp.get('scaffold', '')}\t"
                f"{bp.get('region_type', 'breakpoint')}\t"
                f"{bp.get('center', '')}\t"
                f"{bp.get('start', bp.get('center', ''))}\t"
                f"{bp.get('end', bp.get('center', ''))}\t"
                f"{bp.get('chrom1', '')}\t"
                f"{bp.get('chrom2', '')}\t"
                f"{bp.get('arm1', '')}\t"
                f"{bp.get('arm2', '')}\t"
                f"{bp.get('iscn_notation', '')}\t"
                f"{bp.get('original_scaffold', '')}\n"
            )
    return output_path


def write_tumor_labels_tsv(labels: list, output_path: str) -> str:
    """Write a tumor labels TSV (2 columns, no header)."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        for label in labels:
            f.write(f"{label['scaffold']}\t{label['label']}\n")
    return output_path


def filter_breakpoints_by_scaffolds(breakpoints: list, scaffolds: list) -> list:
    """Filter breakpoints to only those matching the given scaffold names."""
    scaffold_set = set(scaffolds)
    return [bp for bp in breakpoints if bp.get("scaffold") in scaffold_set]
