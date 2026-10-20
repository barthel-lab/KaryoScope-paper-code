# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/__main__.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""CLI entry point for KaryoScope-conductor."""

import argparse
import importlib
import inspect
import pkgutil
import sys
from pathlib import Path

from karyoscope_conductor.pipeline import run_pipeline, load_step


def _discover_steps() -> list:
    """Walk the steps/ package and return all leaf step module paths (e.g. 'iscn.merge_beds')."""
    import karyoscope_conductor.steps as steps_pkg
    steps_dir = Path(steps_pkg.__file__).parent
    discovered = []
    for py_file in sorted(steps_dir.rglob("*.py")):
        if py_file.name == "__init__.py":
            continue
        rel = py_file.relative_to(steps_dir).with_suffix("")
        discovered.append(".".join(rel.parts))
    return discovered


def cmd_list_steps():
    """Print all discoverable steps grouped by domain."""
    steps = _discover_steps()
    if not steps:
        print("No steps found.")
        return
    by_domain = {}
    for s in steps:
        domain, basename = s.split(".", 1) if "." in s else ("(root)", s)
        by_domain.setdefault(domain, []).append(basename)
    for domain in sorted(by_domain):
        print(f"\n{domain}:")
        for basename in sorted(by_domain[domain]):
            print(f"  {domain}.{basename}")


def cmd_describe(step_path: str):
    """Print the docstring and metadata for a step module."""
    try:
        mod = load_step(step_path)
    except ModuleNotFoundError as e:
        print(f"Step '{step_path}' not found: {e}", file=sys.stderr)
        sys.exit(1)
    doc = inspect.getdoc(mod) or "(no docstring)"
    print(f"Step: {step_path}")
    print(f"Module: {mod.__file__}")
    script_name = getattr(mod, "SCRIPT_NAME", None)
    if script_name:
        print(f"Script: {script_name}")
    else:
        basename = step_path.rsplit(".", 1)[-1]
        print(f"Script: KaryoScope_{basename}.py (default)")
    print()
    print(doc)


def main():
    parser = argparse.ArgumentParser(
        prog="karyoscope",
        description="Lightweight YAML-driven workflow runner for KaryoScope analyses",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--config", help="Path to YAML configuration file")
    parser.add_argument("--start-from", dest="start_from", default=None,
                        help="Start pipeline from this step (must be in pipeline:)")
    parser.add_argument("--stop-after", dest="stop_after", default=None,
                        help="Stop pipeline after this step")
    parser.add_argument("--only", nargs="*", default=None,
                        help="Only run analyses whose name/output contains this pattern")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--list-steps", dest="list_steps", action="store_true",
                        help="List all available step modules and exit")
    parser.add_argument("--describe", metavar="STEP", default=None,
                        help="Print info about a step (e.g. --describe iscn.merge_beds)")
    args = parser.parse_args()

    if args.list_steps:
        cmd_list_steps()
        return
    if args.describe:
        cmd_describe(args.describe)
        return
    if not args.config:
        parser.error("--config is required (or use --list-steps / --describe)")

    run_pipeline(args.config, args.start_from, args.stop_after, args.only, args.dry_run)


if __name__ == "__main__":
    main()
