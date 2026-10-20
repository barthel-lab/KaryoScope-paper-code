# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/steps/iscn/assembly_figure.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Karyotype grid figure (TACCAT layout).

Reads:
  - samples: uses sex and assembly_name for the --sample and --sex CLI flags
  - assembly_figure_analyses: list with name, sample, output, figure_args
  - defaults.assembly_figure: passthrough args merged with each analysis's
    figure_args (analysis wins on conflict)

Invokes: KaryoScope_assembly_figure.py (from scripts.iscn).
The script can be overridden per-step via `scripts: { iscn.assembly_figure: <path> }`.
"""

import os

from karyoscope_conductor.common import (
    args_to_cli_flags, run_step, should_run_analysis,
)
from karyoscope_conductor.config import (
    parse_samples, resolve_script, get_merged_args,
)


def run(config: dict, dry_run: bool, only: list | None) -> None:
    analyses = config.get("assembly_figure_analyses", []) or []
    if not analyses:
        return

    samples = parse_samples(config)
    script = resolve_script(__file__, config)

    for analysis in analyses:
        name = analysis["name"]
        output = analysis["output"]
        if not should_run_analysis(name, output, only):
            continue

        print(f"\n  [{name}]")
        sample = samples[analysis["sample"]]
        figure_args = get_merged_args(
            config, analysis.get("figure_args", {}), "assembly_figure"
        )

        cmd = ["python", script, "--output", output]
        if sample.sex:
            cmd.extend(["--sex", sample.sex])
        if sample.assembly_name:
            cmd.extend(["--sample", sample.assembly_name])
        cmd.extend(args_to_cli_flags(figure_args))

        os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
        run_step(cmd, f"assembly_figure {name}", dry_run)
