# Vendored from barthel-lab/KaryoScope-conductor at commit 0732e96adfb4fdfcdc5b2028c10c740c72273ed5, file karyoscope_conductor/config.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""Config loading and shared helpers.

Configs are YAML documents with a `pipeline:` key listing package-style step
paths (e.g. `iscn.merge_beds`). The raw dict is passed to each step module;
steps parse the sections they need via the helpers below.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import yaml

from karyoscope_conductor.common import resolve_path


@dataclass
class Sample:
    """Universal sample representation. Superset of ISCN + BIR fields.

    Steps ignore fields they don't use. ISCN steps read assembly_name/bed_dir/
    sex/role. BIR steps read group/color.
    """
    name: str
    assembly_name: Optional[str] = None
    bed_dir: Optional[str] = None
    sex: Optional[str] = None
    role: Optional[str] = None
    group: Optional[str] = None
    color: Optional[str] = None


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge `override` into `base`. Dicts are merged, scalars/lists replaced."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_raw_config(path: str) -> dict:
    """Load a YAML config and return the raw dict.

    Supports an optional `include:` key (list of paths to merge in). Included
    files are loaded first; the current file's keys recursively override them
    on collision. Paths in `include:` are resolved relative to the current
    file's directory.
    """
    path = resolve_path(path)
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    # Handle include directive (deep merge so nested dicts like `defaults`
    # and `samples` are combined rather than overwritten).
    includes = raw.pop("include", [])
    if includes:
        base_dir = os.path.dirname(path)
        merged = {}
        for inc_path in includes:
            inc_path = inc_path if os.path.isabs(os.path.expanduser(inc_path)) else os.path.join(base_dir, inc_path)
            inc_raw = load_raw_config(inc_path)
            merged = _deep_merge(merged, inc_raw)
        raw = _deep_merge(merged, raw)

    # Validate required top-level keys
    for key in ("name", "pipeline"):
        if key not in raw:
            raise ValueError(f"Config '{path}' missing required key: {key}")

    return raw


def parse_samples(config: dict) -> Dict[str, Sample]:
    """Parse the universal samples section."""
    samples = {}
    for name, data in config.get("samples", {}).items():
        name = str(name)
        data = data or {}
        samples[name] = Sample(
            name=name,
            assembly_name=data.get("assembly_name"),
            bed_dir=data.get("bed_dir"),
            sex=data.get("sex"),
            role=data.get("role"),
            group=data.get("group"),
            color=data.get("color"),
        )
    return samples


def get_scripts_path(config: dict, domain: str) -> str:
    """Resolve the scripts directory for a domain (iscn, analysis, ...)."""
    scripts = config.get("scripts", {}) or {}
    if domain not in scripts:
        raise ValueError(
            f"No scripts directory configured for domain '{domain}'. "
            f"Add to config under 'scripts:' key."
        )
    return resolve_path(scripts[domain])


def resolve_script(module_file: str, config: dict, script_name: Optional[str] = None) -> str:
    """Given __file__ of a step module, resolve the script path to invoke.

    Resolution rules:
      1. Derive domain and basename from the module path:
         .../steps/iscn/merge_beds.py → domain='iscn', basename='merge_beds'
      2. If scripts[f"{domain}.{basename}"] is set, use it as the exact path.
      3. Otherwise use {scripts[domain]}/{script_name or f"KaryoScope_{basename}.py"}
    """
    parts = Path(module_file).parts
    # Expect ".../steps/<domain>/<basename>.py"
    try:
        steps_idx = len(parts) - 1 - parts[::-1].index("steps")
    except ValueError:
        raise ValueError(f"Cannot derive domain/basename from {module_file}")
    domain = parts[steps_idx + 1]
    basename = Path(parts[steps_idx + 2]).stem

    scripts = config.get("scripts", {}) or {}
    override_key = f"{domain}.{basename}"
    if override_key in scripts:
        return resolve_path(scripts[override_key])

    script_filename = script_name or f"KaryoScope_{basename}.py"
    scripts_dir = get_scripts_path(config, domain)
    return os.path.join(scripts_dir, script_filename)


def get_defaults(config: dict, step_basename: str) -> dict:
    """Return the defaults dict for a step basename, handling absent/null sections."""
    return (config.get("defaults", {}) or {}).get(step_basename, {}) or {}


def get_merged_args(config: dict, analysis_args: dict, step_basename: str) -> dict:
    """Merge analysis-level args over defaults for a step basename."""
    defaults = dict(get_defaults(config, step_basename))
    defaults.update(analysis_args or {})
    return defaults


def get_output_prefix(config: dict, analysis: dict) -> str:
    """Return {output_dir-or-output_base_dir}/{analysis.output_prefix}.

    Used by BIR clustering-family steps (cluster_analysis / combine /
    cluster_annotate / cluster_plot) that share a common prefix across steps.
    """
    base = config.get("output_dir") or output_base_dir(config)
    return os.path.join(base, analysis["output_prefix"])


def get_bed_paths_for_samples(config: dict, sample_names: list) -> list:
    """Build BED file paths for a list of sample names using bed_template.

    Uses config["bed_dir"] (global) as the {bed_dir} placeholder, and the
    plain sample name as {sample}. This is the BIR convention where all
    sample BED files live under one shared directory.
    """
    featureset = config.get("featureset", "")
    return [
        expand_bed_template(
            config,
            bed_dir=config.get("bed_dir", ""),
            sample=s,
            featureset=featureset,
        )
        for s in sample_names
    ]


def expand_bed_template(config: dict, **kwargs) -> str:
    """Expand the bed_template with provided keyword arguments."""
    template = config.get("bed_template", "")
    return template.format(
        database=config.get("database", ""),
        **kwargs,
    )


def output_base_dir(config: dict) -> str:
    """Return the output_base_dir, defaulting to 'results'."""
    return config.get("output_base_dir", "results")


def tmp_dir(config: dict) -> str:
    """Return the directory where temporary TSVs are written."""
    path = os.path.join(output_base_dir(config), "tmp")
    os.makedirs(path, exist_ok=True)
    return path
