# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/common.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Shared utilities for step modules."""

import os
import subprocess
import sys
from pathlib import Path


def args_to_cli_flags(args_dict: dict) -> list:
    """Convert a dict of arguments to CLI flag list.

    - key "matrix_type" -> "--matrix-type"
    - bool True -> bare flag (--no-refine)
    - bool False -> omitted
    - None -> omitted (caller used null to explicitly disable a default)
    - list -> flag repeated once per item (--featuresets a --featuresets b)
    - other values -> --key value
    """
    flags = []
    for key, value in args_dict.items():
        flag = "--" + key.replace("_", "-")
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                flags.append(flag)
        elif isinstance(value, list):
            for item in value:
                flags.append(flag)
                flags.append(str(item))
        else:
            flags.append(flag)
            flags.append(str(value))
    return flags


def run_step(cmd: list, step_name: str, dry_run: bool = False):
    """Run a subprocess command with logging."""
    cmd_str = " \\\n    ".join(cmd)
    print(f"\n  $ {cmd_str}")

    if dry_run:
        print("  [dry-run: skipped]")
        return

    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"\n  ERROR: {step_name} failed with exit code {result.returncode}")
        sys.exit(result.returncode)


def should_run_analysis(name: str, output: str, only_patterns: list) -> bool:
    """Check if an analysis should run given --only filter patterns."""
    if not only_patterns:
        return True
    return any(
        pattern.lower() in name.lower() or pattern.lower() in output.lower()
        for pattern in only_patterns
    )


def resolve_path(path: str) -> str:
    """Expand ~ and return an absolute-compatible path string."""
    return str(Path(path).expanduser())
