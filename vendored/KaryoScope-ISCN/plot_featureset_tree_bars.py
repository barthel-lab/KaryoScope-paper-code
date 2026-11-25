#!/usr/bin/env python3
# Vendored from barthel-lab/KaryoScope-ISCN at commit 68a9816baab659c5b710f265c3713844cc2729d2, file plot_featureset_tree_bars.py.
# Copied unmodified except: KS_RESOURCES can be set from the environment (originally fixed to a directory under $HOME);
# draw_panel takes an optional style, and tick_exponents / min_bar_w_for_ticks are added (2026-11-24, Fig. 1b);
# the style can leave out the root label and add an axis title (2026-11-24, round 7), and print display names
# for the node labels (label_name, multigroup_suffix; 2026-11-25). Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""
Visualize KaryoScope featureset hierarchies as vertical trees with
horizontal log-scale bars at each tip, using drawsvg.

Each of the 6 KaryoScope featuresets (chromosome, region, repeat, gene,
acrocentric, subtelomeric) is drawn as a tree showing the priority
hierarchy, with colored bars indicating k-mer counts per feature.

Usage:
    python plot_featureset_tree_bars.py
"""

import argparse
import os
import math
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import drawsvg as draw

# =============================================================================
# Paths
# =============================================================================

PROJECT_DIR = Path(__file__).resolve().parent
KS_RESOURCES = Path(os.environ["KS_RESOURCES"]) if "KS_RESOURCES" in os.environ else Path.home() / "Documents" / "software" / "KaryoScope" / "resources"
HIERARCHY_DIR = KS_RESOURCES / "CHM13"
COLORS_DIR = KS_RESOURCES / "databases" / "KS_human_CHM13"
OUTPUT_DIR = HIERARCHY_DIR

# V2 paths
V2_DB_DIR = KS_RESOURCES / "databases" / "KS_human_CHM13_v2"
V2_OUTPUT_DIR = PROJECT_DIR / "featuresets-hierarchy-v2"

DEFAULT_COLOR = "#CCCCCC"

DATABASE = "KS_human_CHM13"
REFERENCE = "CHM13"

# Abbreviations used in priority filenames ({REFERENCE}_{abbrev}.priority.txt)
FEATURESET_ABBREVS = {
    "chromosome": "chr",
    "region": "reg",
    "repeat": "rep",
    "gene": "gen",
    "acrocentric": "acro",
    "subtelomeric": "subtelo",
}

FEATURESETS = list(FEATURESET_ABBREVS.keys())


def featureset_files(fs_name, database=DATABASE, reference=REFERENCE):
    """Derive filenames for a featureset from database/reference names."""
    abbrev = FEATURESET_ABBREVS[fs_name]
    return {
        "priority": f"{reference}_{abbrev}.priority.txt",
        "dbinfo": f"{database}.{fs_name}.database.info.txt",
        "colors": f"{database}.{fs_name}.colors.txt",
    }

# Color overrides for features missing or needing different colors
COLOR_OVERRIDES = {
    "rDNA_specific": "#964B00",          # brown (like noncentromeric)
    "srpRNA_specific": "#90D5FF",        # like other RNA
    "mixedAlpha_specific": "#fa0000",    # like HOR (active)
    "Satellite_specific": "#008B8B",     # teal
    "Simple_repeat_specific": "#DAA520", # gold
    "Low_complexity_specific": "#6B8E23",# olive
    "Noninterspersed_multigroup1": "#B0C4DE",  # match other multigroup nodes
}

# Layout constants
FONT_SIZE = 14
FONT_FAMILY = "Helvetica"
ROW_H = 16          # row height in px
BAR_H = 7           # bar thickness — matches font cap height
MAX_BAR_W = 60      # max bar width in px
INDENT = 14         # x indent per depth level
TICK_LEN = 8        # horizontal connector tick
TEXT_GAP = 6        # gap between text end and bar start
CHAR_W = 8.0        # approx character width at font_size=14
LINE_COLOR = "#000000"
LINE_W = 1


# =============================================================================
# Data model
# =============================================================================


@dataclass
class TreeNode:
    name: str
    priority: int = 0
    parent_name: str | None = None
    kmer_count: int = 0
    color: str = DEFAULT_COLOR
    children: list["TreeNode"] = field(default_factory=list)
    is_leaf: bool = True


# =============================================================================
# Parsing
# =============================================================================


def parse_priority_file(filepath: Path) -> list[tuple[str, int, str]]:
    """Parse priority.txt -> [(child, priority, parent), ...]"""
    entries = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                child = parts[0]
                priority = int(parts[1])
                parent = parts[2]
                entries.append((child, priority, parent))
    return entries


def parse_database_info(filepath: Path) -> dict[str, int]:
    """Parse database.info -> {feature_suffix: count}"""
    counts = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                name = parts[0]
                if name == "TOTAL":
                    continue
                count = int(parts[1].replace(",", ""))
                counts[name] = count
    return counts


def parse_color_file(filepath: Path) -> dict[str, str]:
    """Parse colors.txt -> {feature_suffix: hex_color}"""
    colors = {}
    with open(filepath) as f:
        for i, line in enumerate(f):
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                if i == 0 and parts[0].lower() == "feature":
                    continue
                colors[parts[0]] = parts[1]
    return colors


def build_tree(
    priority_tuples: list[tuple[str, int, str]],
    kmer_counts: dict[str, int],
    colors: dict[str, str],
) -> TreeNode:
    """Build tree from parsed data, returning the root 'categorized' node."""
    parent_names = {parent for _, _, parent in priority_tuples}

    nodes: dict[str, TreeNode] = {}

    root = TreeNode(name="categorized")
    root.is_leaf = "categorized" in parent_names
    nodes["categorized"] = root

    for child, priority, parent in priority_tuples:
        node = TreeNode(name=child, priority=priority, parent_name=parent)
        node.is_leaf = child not in parent_names
        nodes[child] = node

    for name, node in nodes.items():
        if node.is_leaf:
            suffix = f"{name}_specific"
        else:
            suffix = f"{name}_multigroup1"
        node.kmer_count = kmer_counts.get(suffix, 0)
        node.color = colors.get(suffix, DEFAULT_COLOR)

    for child, priority, parent in priority_tuples:
        if parent in nodes:
            nodes[parent].children.append(nodes[child])
            nodes[parent].is_leaf = False

    for node in nodes.values():
        node.children.sort(key=lambda n: n.priority)

    return root


def load_all_featuresets() -> dict[str, TreeNode]:
    """Load all 6 v1 featuresets, returning {name: root_node}."""
    featuresets = {}
    for fs_name in FEATURESETS:
        files = featureset_files(fs_name)
        priority_tuples = parse_priority_file(HIERARCHY_DIR / files["priority"])
        kmer_counts = parse_database_info(COLORS_DIR / files["dbinfo"])
        colors = parse_color_file(COLORS_DIR / files["colors"])
        featuresets[fs_name] = build_tree(priority_tuples, kmer_counts, colors)
    return featuresets


# =============================================================================
# V2 support
# =============================================================================


def parse_v2_priority_file(filepath: Path, featureset: str) -> list[tuple[str, int, str]]:
    """Parse v2 combined priority file, filtering by featureset.
    Format: feature_set<tab>child<tab>priority<tab>parent"""
    entries = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                if parts[0].lower() == "feature_set":
                    continue
                if parts[0] == featureset:
                    child, priority, parent = parts[1], int(parts[2]), parts[3]
                    entries.append((child, priority, parent))
    return entries


def build_tree_v2(
    priority_tuples: list[tuple[str, int, str]],
    kmer_counts: dict[str, int],
    colors: dict[str, str],
) -> TreeNode:
    """Build tree from v2 data (bare feature names, no _specific/_multigroup1 suffixes)."""
    parent_names = {parent for _, _, parent in priority_tuples}

    nodes: dict[str, TreeNode] = {}

    root = TreeNode(name="categorized")
    root.is_leaf = "categorized" in parent_names
    nodes["categorized"] = root

    for child, priority, parent in priority_tuples:
        node = TreeNode(name=child, priority=priority, parent_name=parent)
        node.is_leaf = child not in parent_names
        nodes[child] = node

    for name, node in nodes.items():
        node.kmer_count = kmer_counts.get(name, 0)
        node.color = colors.get(name, DEFAULT_COLOR)

    for child, priority, parent in priority_tuples:
        if parent in nodes:
            nodes[parent].children.append(nodes[child])
            nodes[parent].is_leaf = False

    for node in nodes.values():
        node.children.sort(key=lambda n: n.priority)

    return root


def load_all_featuresets_v2() -> dict[str, TreeNode]:
    """Load v2 featuresets from combined priority file."""
    priority_file = V2_DB_DIR / "KS_human_CHM13_v2.priority.txt"
    all_fs_names = set()
    with open(priority_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[0].lower() != "feature_set":
                all_fs_names.add(parts[0])

    # V2 dbinfo file naming: KS_human_CHM13.CHM13.{featureset}.database.info.txt
    featuresets = {}
    for fs_name in all_fs_names:
        priority_tuples = parse_v2_priority_file(priority_file, fs_name)

        # Load k-mer counts from v2 database.info file (bare feature names)
        dbinfo_path = V2_DB_DIR / f"KS_human_CHM13.{fs_name}.database.info.txt"
        if dbinfo_path.exists():
            kmer_counts = parse_database_info(dbinfo_path)
        else:
            # Fallback to 1M placeholder
            all_names = {child for child, _, _ in priority_tuples}
            all_names.add("categorized")
            kmer_counts = {name: 1_000_000 for name in all_names}

        # Load colors
        colors_path = V2_DB_DIR / f"KS_human_CHM13_v2.{fs_name}.colors.txt"
        colors = parse_color_file(colors_path) if colors_path.exists() else {}
        featuresets[fs_name] = build_tree_v2(priority_tuples, kmer_counts, colors)
    return featuresets


# =============================================================================
# Formatting helpers
# =============================================================================


def format_count(n: int) -> str:
    """Human-readable k-mer count."""
    if n >= 1_000_000_000:
        return f"{n / 1_000_000_000:.1f}B"
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n / 1_000:.1f}K"
    return str(n)


# =============================================================================
# Flatten tree
# =============================================================================


def _flatten(node: TreeNode, depth: int = 0) -> list[dict]:
    """DFS pre-order flatten, recording depth."""
    rows = [{"name": node.name, "depth": depth, "node": node,
             "is_leaf": not node.children}]
    for child in node.children:
        rows.extend(_flatten(child, depth + 1))
    return rows


# =============================================================================
# Draw a single featureset panel, returns (group, width, height)
# =============================================================================


def draw_panel(root: TreeNode, featureset: str, style: dict | None = None) -> tuple:
    """Draw one featureset tree+bars, returning (draw.Group, width, height).

    style (KaryoScope-paper-code, 2026-11-24) overrides the layout for a figure at final size:
    font_size, root_font, count_font, tick_font, row_h, bar_h, max_bar_w, indent, tick_len,
    text_gap, line_w, font_family, multigroup (False drops the " (multigroup)" suffix) and
    text_width (a function (text, font_size, bold) -> width, replacing the CHAR_W estimate),
    root_label (False leaves out the bold root label, the tree's header; the layout is unchanged),
    label_name (a function name -> printed label, replacing name.replace("_", " ")),
    multigroup_suffix (the text after an inner node's count, default " (multigroup)"),
    parent_bar_frac (inner nodes' bars drawn this fraction of BAR_H tall, centred; default 1),
    axis_title (a list of (text, italic) parts, set on the tick-label baseline and ending just left
    of the first tick label), axis_title_font and axis_title_gap (default half the tick font).
    Without it the drawing is unchanged."""
    st = style or {}
    FONT_SIZE = st.get("font_size", globals()["FONT_SIZE"])
    ROOT_FONT = st.get("root_font", FONT_SIZE + 2)
    COUNT_FONT = st.get("count_font", FONT_SIZE - 1.5)
    TICK_FONT = st.get("tick_font", 6.5)
    ROW_H = st.get("row_h", globals()["ROW_H"])
    BAR_H = st.get("bar_h", globals()["BAR_H"])
    MAX_BAR_W = st.get("max_bar_w", globals()["MAX_BAR_W"])
    INDENT = st.get("indent", globals()["INDENT"])
    TICK_LEN = st.get("tick_len", globals()["TICK_LEN"])
    TEXT_GAP = st.get("text_gap", globals()["TEXT_GAP"])
    LINE_W = st.get("line_w", globals()["LINE_W"])
    FONT_FAMILY = st.get("font_family", globals()["FONT_FAMILY"])
    text_width = st.get("text_width")
    multigroup = st.get("multigroup", True)
    MIN_W = st.get("min_w", 200)
    root_label = st.get("root_label", True)
    label_name = st.get("label_name")
    multigroup_suffix = st.get("multigroup_suffix", " (multigroup)")
    parent_bar_frac = st.get("parent_bar_frac", 1.0)
    axis_title = st.get("axis_title")

    rows = _flatten(root)
    # Rename root "categorized" to the featureset name
    for r in rows:
        if r["name"] == "categorized":
            r["name"] = featureset

    n = len(rows)
    if n == 0:
        return draw.Group(), 100, 30

    # Compute log scale
    all_counts = [r["node"].kmer_count for r in rows if r["node"].kmer_count > 0]
    log_max = math.log10(max(all_counts)) if all_counts else 1

    # Assign y positions
    title_h = 4
    for i, r in enumerate(rows):
        r["y"] = title_h + i * ROW_H + ROW_H / 2

    # Build parent map (row index -> parent row index)
    parent_idx = {}
    for i, r in enumerate(rows):
        if r["depth"] == 0:
            parent_idx[i] = None
        else:
            for j in range(i - 1, -1, -1):
                if rows[j]["depth"] == r["depth"] - 1:
                    parent_idx[i] = j
                    break

    children_of: dict[int, list[int]] = {i: [] for i in range(n)}
    for i, pi in parent_idx.items():
        if pi is not None:
            children_of[pi].append(i)

    g = draw.Group()

    # Draw tree connectors
    for pi, cis in children_of.items():
        if not cis:
            continue
        parent_depth = rows[pi]["depth"]
        parent_y = rows[pi]["y"]

        vx = (parent_depth + 1) * INDENT

        last_child_y = rows[cis[-1]]["y"]

        start_y = parent_y + ROW_H * 0.5
        g.append(draw.Line(vx, start_y, vx, last_child_y,
                           stroke=LINE_COLOR, stroke_width=LINE_W))

        for ci in cis:
            cy = rows[ci]["y"]
            g.append(draw.Line(vx, cy, vx + TICK_LEN, cy,
                               stroke=LINE_COLOR, stroke_width=LINE_W))

    # --- First pass: compute text positions and find aligned bar_left ---
    text_info = []
    for i, r in enumerate(rows):
        depth = r["depth"]
        name = r["name"]
        display_name = label_name(name) if label_name and not (depth == 0 and i == 0) else name.replace("_", " ")
        r["display_name"] = display_name
        is_root = (depth == 0 and i == 0)
        if depth == 0:
            tx = 0
        else:
            tx = depth * INDENT + TICK_LEN + 3
        if text_width:
            text_end = tx + text_width(display_name, ROOT_FONT if is_root else FONT_SIZE, is_root)
        else:
            char_w = CHAR_W * (1.15 if is_root else 1.0)
            text_end = tx + len(display_name) * char_w
        text_info.append({"tx": tx, "text_end": text_end, "is_root": is_root})

    bar_left = max(t["text_end"] for t in text_info) + TEXT_GAP + 4
    dot_gap = 3

    # --- Draw axis at the root row level ---
    root_y = rows[0]["y"]
    axis_y = root_y

    max_exp = int(math.floor(log_max))
    if max_exp <= 4:
        step = 1
    elif max_exp <= 7:
        step = 2
    else:
        step = 3
    tick_exps = list(range(0, max_exp + 1, step))
    if max_exp not in tick_exps:
        if max_exp - tick_exps[-1] >= step * 0.6:
            tick_exps.append(max_exp)

    g.append(draw.Line(bar_left, axis_y, bar_left + MAX_BAR_W, axis_y,
                       stroke="#000000", stroke_width=0.7))
    for exp in tick_exps:
        tx_pos = bar_left + (exp / log_max) * MAX_BAR_W if log_max > 0 else bar_left
        g.append(draw.Line(tx_pos, axis_y - 3, tx_pos, axis_y,
                           stroke="#000000", stroke_width=0.7))
        tick_val = 10 ** exp
        tick_str = format_count(tick_val)
        g.append(draw.Text(
            tick_str, TICK_FONT,
            x=tx_pos, y=axis_y - 5,
            fill="#000000", font_family=FONT_FAMILY,
            text_anchor="middle",
        ))
    if axis_title:
        # what the axis counts, on the tick-label row, ending just left of the first tick label
        tfont = st.get("axis_title_font", TICK_FONT)
        gap = st.get("axis_title_gap", 0.5 * TICK_FONT)
        width = (lambda t, f: text_width(t, f, False)) if text_width else (lambda t, f: len(t) * CHAR_W * f / FONT_SIZE)
        first = format_count(10 ** tick_exps[0])
        x_end = bar_left - width(first, TICK_FONT) / 2 - gap
        title_w = sum(width(t, tfont) for t, _ in axis_title)
        if x_end - title_w < 0:
            raise ValueError(f"axis title ({title_w:.1f}) does not fit left of the first tick ({x_end:.1f})")
        # one left-anchored text per part at its own x (renderers differ on text-anchor="end" with an
        # italic <tspan>: cairosvg puts the tspans after the end point)
        x = x_end - title_w
        for t, italic in axis_title:
            g.append(draw.Text(t, tfont, x=x, y=axis_y - 5, fill="#000000", font_family=FONT_FAMILY,
                               font_style="italic" if italic else "normal", text_anchor="start"))
            x += width(t, tfont)

    # --- Second pass: draw labels, dotted leaders, bars, counts ---
    max_right = 0
    for i, r in enumerate(rows):
        y = r["y"]
        node = r["node"]
        count = node.kmer_count
        color = node.color
        is_leaf = r["is_leaf"]
        ti = text_info[i]
        tx = ti["tx"]
        text_end = ti["text_end"]
        is_root = ti["is_root"]

        display_name = r["display_name"]
        if root_label or not is_root:
            g.append(draw.Text(
                display_name, ROOT_FONT if is_root else FONT_SIZE,
                x=tx, y=y + FONT_SIZE * 0.35,
                fill="#000000" if is_root else "#333333",
                font_family=FONT_FAMILY,
                font_weight="bold" if is_root else "normal",
            ))

        if is_root:
            max_right = max(max_right, bar_left + MAX_BAR_W + 10)
            continue

        leader_start = text_end + dot_gap
        leader_end = bar_left - 2
        if leader_end > leader_start + 2:
            g.append(draw.Line(
                leader_start, y, leader_end, y,
                stroke="#999999", stroke_width=0.5,
                stroke_dasharray="1.5,2",
            ))

        if count > 0:
            log_val = math.log10(count)
            bw = max(2, (log_val / log_max) * MAX_BAR_W) if log_max > 0 else 5
        else:
            bw = 1.5

        bh = BAR_H * parent_bar_frac if not is_leaf and parent_bar_frac != 1.0 else BAR_H
        g.append(draw.Rectangle(
            bar_left, y - bh / 2, bw, bh,
            fill=color, stroke="none",
        ))

        count_str = format_count(count) if count > 0 else "0"
        suffix = multigroup_suffix if not is_leaf and count > 0 and multigroup else ""
        g.append(draw.Text(
            count_str + suffix, COUNT_FONT,
            x=bar_left + bw + 4, y=y + COUNT_FONT * 0.35,
            fill="#555555", font_family=FONT_FAMILY,
        ))

        if text_width:
            right_edge = bar_left + bw + 4 + text_width(count_str + suffix, COUNT_FONT, False)
        else:
            right_edge = bar_left + bw + 4 + len(count_str + suffix) * 4.5 + 10
        if right_edge > max_right:
            max_right = right_edge

    total_h = title_h + n * ROW_H + 4
    total_w = max(max_right, MIN_W)

    return g, total_w, total_h


def tick_exponents(root: TreeNode) -> tuple:
    """The axis ticks draw_panel draws for this tree: (log_max, [exponents])."""
    counts = [r["node"].kmer_count for r in _flatten(root) if r["node"].kmer_count > 0]
    log_max = math.log10(max(counts)) if counts else 1
    max_exp = int(math.floor(log_max))
    step = 1 if max_exp <= 4 else 2 if max_exp <= 7 else 3
    exps = list(range(0, max_exp + 1, step))
    if max_exp not in exps and max_exp - exps[-1] >= step * 0.6:
        exps.append(max_exp)
    return log_max, exps


def min_bar_w_for_ticks(root: TreeNode, tick_font: float, text_width, gap: float) -> float:
    """Shortest MAX_BAR_W (tree units) at which draw_panel's centred tick labels, at tick_font,
    are at least `gap` apart (KaryoScope-paper-code, 2026-11-24)."""
    log_max, exps = tick_exponents(root)
    need = 0.0
    for e0, e1 in zip(exps, exps[1:]):
        w0 = text_width(format_count(10 ** e0), tick_font, False)
        w1 = text_width(format_count(10 ** e1), tick_font, False)
        need = max(need, (w0 / 2 + w1 / 2 + gap) * log_max / (e1 - e0))
    return need


# =============================================================================
# Main
# =============================================================================


def _apply_color_overrides(node: TreeNode) -> None:
    """Apply color overrides to tree nodes."""
    suffix = f"{node.name}_specific" if node.is_leaf else f"{node.name}_multigroup1"
    if suffix in COLOR_OVERRIDES:
        node.color = COLOR_OVERRIDES[suffix]
    for child in node.children:
        _apply_color_overrides(child)


def _save_with_png(d, out_svg):
    """Save SVG and attempt PNG conversion."""
    d.save_svg(str(out_svg))
    print(f"Saved: {out_svg}")
    out_png = out_svg.with_suffix(".png")
    from karyoplot.svg.export import svg_to_png
    if svg_to_png(str(out_svg), str(out_png), dpi=300, quiet=True):
        print(f"Saved: {out_png}")


def _draw_individual(featuresets, fs_name):
    """Draw a single featureset as a standalone SVG."""
    grp, pw, ph = draw_panel(featuresets[fs_name], fs_name)
    margin = 10
    d = draw.Drawing(pw + margin * 2, ph + margin * 2, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, pw + margin * 2, ph + margin * 2, fill="white"))
    container = draw.Group(transform=f"translate({margin},{margin})")
    container.append(grp)
    d.append(container)
    out_svg = OUTPUT_DIR / f"{fs_name}_tree.svg"
    _save_with_png(d, out_svg)


def _draw_grid(featuresets):
    """Draw the combined 2x3 grid of all featuresets."""
    order = ["chromosome", "region", "repeat", "gene", "acrocentric", "subtelomeric"]
    panels = {fs: draw_panel(featuresets[fs], fs) for fs in order}

    col_gap = 30
    row_gap = 20
    row1, row2 = order[:3], order[3:]

    col_widths = [max(panels[row1[c]][1], panels[row2[c]][1]) for c in range(3)]
    row_heights = [
        max(panels[fs][2] for fs in row1),
        max(panels[fs][2] for fs in row2),
    ]

    total_w = sum(col_widths) + col_gap * 2 + 20
    total_h = sum(row_heights) + row_gap + 40

    d = draw.Drawing(total_w, total_h, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, total_w, total_h, fill="white"))
    d.append(draw.Text(
        "KaryoScope Featureset Hierarchies", 14,
        x=total_w / 2, y=16,
        fill="#000000", font_family=FONT_FAMILY, font_weight="bold",
        text_anchor="middle",
    ))

    y_offset = 36
    for row_idx, fs_row in enumerate([row1, row2]):
        x_offset = 10
        for col_idx, fs_name in enumerate(fs_row):
            grp, pw, ph = panels[fs_name]
            container = draw.Group(transform=f"translate({x_offset},{y_offset})")
            container.append(grp)
            d.append(container)
            x_offset += col_widths[col_idx] + col_gap
        y_offset += row_heights[row_idx] + row_gap

    _save_with_png(d, OUTPUT_DIR / "featureset_tree_bars.svg")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize KaryoScope featureset hierarchies as trees with bars")
    parser.add_argument("--featuresets", default=None,
                        help="Comma-separated featuresets to plot individually "
                             "(default: all in a combined grid)")
    parser.add_argument("--version", default="v1", choices=["v1", "v2"],
                        help="Database version (default: v1)")
    args = parser.parse_args()

    if args.version == "v2":
        all_fs = load_all_featuresets_v2()
        # Use v2 output directory
        global OUTPUT_DIR
        V2_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        OUTPUT_DIR = V2_OUTPUT_DIR
    else:
        all_fs = load_all_featuresets()

    for root in all_fs.values():
        _apply_color_overrides(root)

    if args.featuresets:
        for fs_name in args.featuresets.split(","):
            fs_name = fs_name.strip()
            if fs_name not in all_fs:
                print(f"Warning: unknown featureset '{fs_name}', skipping")
                continue
            _draw_individual(all_fs, fs_name)
    else:
        _draw_grid(all_fs)


if __name__ == "__main__":
    main()
