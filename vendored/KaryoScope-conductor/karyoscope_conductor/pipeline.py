# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/pipeline.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Pipeline orchestrator.

Loads a YAML config, reads its `pipeline:` list of step paths (dotted notation),
and dispatches each step to the corresponding module under
`karyoscope_conductor.steps.<path>`.
"""

import importlib
from typing import Optional, List

from karyoscope_conductor.config import load_raw_config


def load_step(step_path: str):
    """Given 'iscn.merge_beds', import karyoscope_conductor.steps.iscn.merge_beds."""
    return importlib.import_module(f"karyoscope_conductor.steps.{step_path}")


def run_pipeline(
    config_path: str,
    start_from: Optional[str] = None,
    stop_after: Optional[str] = None,
    only: Optional[List[str]] = None,
    dry_run: bool = False,
):
    config = load_raw_config(config_path)
    pipeline = config.get("pipeline") or []

    if not pipeline:
        print("Warning: empty pipeline in config.")
        return

    if start_from and start_from not in pipeline:
        raise ValueError(f"--start-from '{start_from}' not in pipeline: {pipeline}")
    if stop_after and stop_after not in pipeline:
        raise ValueError(f"--stop-after '{stop_after}' not in pipeline: {pipeline}")

    start_idx = pipeline.index(start_from) if start_from else 0
    stop_idx = pipeline.index(stop_after) if stop_after else len(pipeline) - 1

    steps_to_run = pipeline[start_idx:stop_idx + 1]

    print("=" * 60)
    print(f"  {config.get('name', '(unnamed)')}")
    print("=" * 60)
    print(f"  Pipeline: {' -> '.join(pipeline)}")
    if start_from or stop_after:
        print(f"  Running: {' -> '.join(steps_to_run)}")
    if only:
        print(f"  Filter: {only}")

    for step_path in steps_to_run:
        step = load_step(step_path)
        print(f"\n{'=' * 60}")
        print(f"  Step: {step_path}")
        print("=" * 60)
        step.run(config, dry_run, only)

    print(f"\n{'=' * 60}")
    print("  Pipeline complete!")
    print("=" * 60)
