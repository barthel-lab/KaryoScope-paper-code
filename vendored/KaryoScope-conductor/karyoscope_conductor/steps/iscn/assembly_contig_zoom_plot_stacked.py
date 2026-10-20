# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/steps/iscn/assembly_contig_zoom_plot_stacked.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Stacked multi-sample zoom plots.

Same underlying script as assembly_contig_zoom_plot, but invoked with
repeated --sample/--bed-dir/--regions triples plus --stacked-output.

Reads:
  - samples
  - breakpoint_groups
  - assembly_contig_zoom_plot_stacked_analyses: list with entries=[{sample, breakpoint_group}]
  - defaults.assembly_contig_zoom_plot_stacked: passthrough args
  - colors_dir, database
"""

import os

from karyoscope_conductor.common import (
    args_to_cli_flags, run_step, should_run_analysis, resolve_path,
)
from karyoscope_conductor.config import (
    parse_samples, resolve_script, get_merged_args, tmp_dir,
)
from karyoscope_conductor.tsv_generator import write_regions_tsv


# Reuses the same script as the non-stacked variant
SCRIPT_NAME = "KaryoScope_assembly_contig_zoom_plot.py"


def run(config: dict, dry_run: bool, only: list | None) -> None:
    analyses = config.get("assembly_contig_zoom_plot_stacked_analyses", []) or []
    if not analyses:
        return

    samples = parse_samples(config)
    breakpoint_groups = config.get("breakpoint_groups", {}) or {}
    script = resolve_script(__file__, config, script_name=SCRIPT_NAME)
    colors_dir = resolve_path(config.get("colors_dir", ""))
    database = config.get("database", "")

    for analysis in analyses:
        name = analysis["name"]
        output = analysis["output"]
        if not should_run_analysis(name, output, only):
            continue

        print(f"\n  [{name}]")
        stacked_args = get_merged_args(
            config, analysis.get("stacked_args", {}), "assembly_contig_zoom_plot_stacked"
        )
        safe_name = name.replace(" ", "_").replace("/", "_")

        cmd = ["python", script]

        for i, entry in enumerate(analysis["entries"]):
            sample = samples[entry["sample"]]
            breakpoints = breakpoint_groups.get(entry["breakpoint_group"], [])
            regions_path = os.path.join(
                tmp_dir(config), f"{safe_name}_entry{i}_{entry['sample']}_regions.tsv"
            )
            if not dry_run:
                write_regions_tsv(breakpoints, regions_path)
            cmd.extend(["--sample", sample.assembly_name])
            cmd.extend(["--bed-dir", sample.bed_dir])
            cmd.extend(["--regions", regions_path])

        cmd.extend([
            "--colors-dir", colors_dir,
            "--database", database,
            "--output-dir", os.path.dirname(output) or ".",
            "--stacked-output", output,
            *args_to_cli_flags(stacked_args),
        ])

        os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
        run_step(cmd, f"stacked {name}", dry_run)
