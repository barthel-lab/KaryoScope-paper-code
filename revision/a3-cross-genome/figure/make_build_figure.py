#!/usr/bin/env python3
"""Assemble the build figure: panels A and B from the hand-authored draft,
panel C generated here from the A3 cross-genome run.

Panel C shows the Ey15-2 transfer beside the Col-CEN self-annotation control,
plus the label composition that the genome-view renderer cannot show.

Why the composition bars are not decoration: the genome view quantises to 4 px
blocks (~138 kb at the native Ey15-2 scale) and its legend carries only the four
leaf labels it drew. Every minority label -- arm, novel, categorized,
centromeric -- is not drawn. A bare side-by-side would therefore
read as a flawless transfer while the Results text reports 10.3% falling back to
a parent node.

    python3 scripts/make_build_figure.py [--a3-dir DIR]

Panel C reads the A3 cross-genome run directory (the two karyotype renders and
results/label_composition.tsv). It is located by --a3-dir, else the KS_A3_DIR
environment variable, else <repo>/a3_crossgenome.

Writes figures/fig7/build_workflow.svg.

Text sizes are SVG user units on the 1400-unit-wide canvas. The figure is
printed 183 mm wide, so 1 unit = 0.3705 pt; the FS_* constants keep body text
at >= 5 pt and panel letters at 8 pt at that size.
"""

import os
import re
import csv
import argparse
import collections
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DRAFT = ROOT / "figures" / "fig7" / "build_workflow_panelsAB.svg"
OUT = ROOT / "figures" / "fig7" / "build_workflow.svg"
A3 = Path(os.environ.get("KS_A3_DIR", ROOT / "a3_crossgenome"))


def set_a3(a3):
    """Point panel C's three inputs at an A3 run directory."""
    global A3, EY_SVG, CC_SVG, COMP
    A3 = Path(a3)
    EY_SVG = A3 / "karyotype_ey15_2" / \
        "ey15-2.HKS_arabidopsis_ColCEN.genome.region.smoothed.karyotype.svg"
    CC_SVG = A3 / "karyotype_colcen" / \
        "Col-CEN_v1.2.HKS_arabidopsis_ColCEN.genome.region.smoothed.karyotype.svg"
    COMP = A3 / "results" / "label_composition.tsv"


set_a3(A3)

# Scale each source render was drawn at, from the A3 handoff. The two differ by
# ~6%, which is exactly why they cannot be placed side by side as rendered.
KBP_PER_PX = {"ey": 34.4, "cc": 32.5}

# ------------------------------------------------------------ font sizes -----
# User units; x 0.3705 = pt at the printed width of 183 mm.
FS_LETTER = 21.6    # 8.0 pt, panel letter
FS_TITLE = 18.5     # 6.9 pt, panel title
FS_GROUP = 16       # 5.9 pt, genome group labels
FS_BODY = 14        # 5.2 pt, labels, legends, headers
FS_SMALL = 13.5     # 5.0 pt, tick labels, bar values, parent-label names
CHAR_W = 8.5 * FS_BODY / 13   # legend advance per character (was 8.5 at 13)

# ---------------------------------------------------------------- canvas -----
W = 1400
H = 1290
PANEL_C_TOP = 684

# Karyotype block (left)
KX0 = 92            # left edge of the chromosome tracks
KW = 700            # width allotted to the longest chromosome
BAR_H = 14
ROW_PITCH = 22

# Composition block (right)
CX0 = 856
CW = 500

# The tool's own palette, so panel C matches KaryoScope output.
LEAF_COLOR = {
    "q_arm": "#732EB8",
    "p_arm": "#D14747",
    "CEN180": "#73B82E",
    "cen_gap": "#47D1D1",
}
# Non-leaf labels get greys so "generic" reads as generic; novel gets amber.
NONLEAF_COLOR = {
    "arm": "#8A8A8A",
    "centromeric": "#ADADAD",
    "categorized": "#CFCFCF",
    "novel": "#E8A33D",
}
LEAF_ORDER = ["q_arm", "p_arm", "CEN180", "cen_gap"]
NONLEAF_ORDER = ["arm", "novel", "categorized", "centromeric"]

# The one generic label per feature set, as fixed by the A3 README. The repeat
# hierarchy is not reconstructed here -- only the label the README names is
# quoted, so nothing depends on guessing which repeat labels are internal nodes.
PARENT_LABEL = {
    "chromosome": "nuclear",
    "region": "arm",
    "gene": "categorized",
    "repeat": "Transposable_Element",
}
SET_ORDER = ["chromosome", "region", "gene", "repeat"]

RECT_RE = re.compile(
    r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" '
    r'height="([\d.]+)" fill="(#[0-9A-Fa-f]{6})"'
)


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def read_karyotype(path):
    """Pull the annotation tracks out of a genome-view karyotype render.

    The tracks are the 6 px-wide rect column per chromosome; everything else in
    the file is the title, the chromosome tick labels, the scale bar and the
    legend. Returns one list of (y, height, fill) per chromosome, left to right.
    """
    text = path.read_text()
    rects = RECT_RE.findall(text)
    track = [r for r in rects if abs(float(r[2]) - 6.0) < 1e-6]
    if not track:
        raise SystemExit(f"no 6px annotation track found in {path}")
    cols = collections.OrderedDict()
    for x, y, _w, h, fill in track:
        cols.setdefault(float(x), []).append((float(y), float(h), fill))
    out = []
    for x in sorted(cols):
        blocks = sorted(cols[x])
        out.append(blocks)
    names = re.findall(r'<text x="[\d.]+" y="58"[^>]*>([^<]+)<', text)
    return out, names


# Tiled genome lengths, from the A3 label-composition README. Used to recover
# true precision for labels whose pct column rounds to zero -- the chromosome
# set's novel fraction is 105 bp, which the table stores as 0.0001.
TILED_BP = {"ey15_2": 136135221, "colcen": 131559526}


def load_composition():
    comp = collections.defaultdict(dict)
    with COMP.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            run, fset = row["run"], row["feature_set"]
            pct = float(row["bp"]) / TILED_BP[run] * 100.0
            comp[(run, fset)][row["label"]] = pct
    return comp


def fmt_pct(v):
    """Two decimals normally, more only where two would read as a flat zero."""
    if v == 0:
        return "0"
    if v >= 0.01:
        return f"{v:.2f}"
    if v >= 0.0001:
        return f"{v:.4f}"
    return f"{v:.5f}"


def karyotype_group(blocks, kbp_per_px, kb_per_px_target, x0, y0, names):
    """Re-project a vertical karyotype as horizontal bars at a shared scale.

    Scaling is done on the source geometry rather than by re-deriving base
    pairs, so the only quantisation is the renderer's own 4 px blocks.
    """
    out = []
    factor = kbp_per_px / kb_per_px_target
    for i, chrom in enumerate(blocks):
        y = y0 + i * ROW_PITCH
        top = chrom[0][0]
        name = names[i] if i < len(names) else f"Chr{i + 1}"
        out.append(
            f'  <text x="{x0 - 8}" y="{y + BAR_H - 3:.0f}" font-size="{FS_BODY}" '
            f'fill="#555555" text-anchor="end">{esc(name)}</text>'
        )
        for by, bh, fill in chrom:
            bx = x0 + (by - top) * factor
            bw = bh * factor
            out.append(
                f'  <rect x="{bx:.2f}" y="{y}" width="{bw:.2f}" '
                f'height="{BAR_H}" fill="{fill}"/>'
            )
        total = sum(b[1] for b in chrom) * factor
        out.append(
            f'  <rect x="{x0}" y="{y}" width="{total:.2f}" height="{BAR_H}" '
            f'fill="none" stroke="#666666" stroke-width="0.6"/>'
        )
    return out


def stacked_bar(pcts, x0, y, width, height):
    """One 100%-stacked composition bar: leaves first, then the non-leaf mass."""
    out = []
    cursor = float(x0)
    for label in LEAF_ORDER + NONLEAF_ORDER:
        pct = pcts.get(label, 0.0)
        if pct <= 0:
            continue
        w = width * pct / 100.0
        color = LEAF_COLOR.get(label) or NONLEAF_COLOR[label]
        out.append(
            f'  <rect x="{cursor:.2f}" y="{y}" width="{w:.3f}" '
            f'height="{height}" fill="{color}"/>'
        )
        cursor += w
    out.append(
        f'  <rect x="{x0}" y="{y}" width="{width}" height="{height}" '
        f'fill="none" stroke="#666666" stroke-width="0.8"/>'
    )
    return out


def nonleaf_total(pcts):
    return sum(pcts.get(k, 0.0) for k in NONLEAF_ORDER)


def build_panel_c(comp):
    ey_blocks, ey_names = read_karyotype(EY_SVG)
    cc_blocks, cc_names = read_karyotype(CC_SVG)

    # Longest chromosome across both genomes sets the shared scale.
    ey_longest = max(sum(b[1] for b in c) for c in ey_blocks) * KBP_PER_PX["ey"]
    cc_longest = max(sum(b[1] for b in c) for c in cc_blocks) * KBP_PER_PX["cc"]
    longest_kb = max(ey_longest, cc_longest)
    kb_per_px = longest_kb / KW

    s = []
    top = PANEL_C_TOP
    s.append(f'  <text x="20" y="{top + 26}" font-size="{FS_LETTER}" font-weight="bold" '
             f'fill="#222222">c</text>')
    s.append(f'  <text x="50" y="{top + 26}" font-size="{FS_TITLE}" fill="#222222">'
             f'Annotating a second accession with the database</text>')

    # ---- left: the two karyotypes, matched scale --------------------------
    ey_y = top + 60
    s.append(f'  <text x="{KX0 - 52}" y="{ey_y - 12}" font-size="{FS_GROUP}" '
             f'font-weight="bold" fill="#333333">Ey15-2</text>')
    s += karyotype_group(ey_blocks, KBP_PER_PX["ey"], kb_per_px, KX0, ey_y,
                         ey_names)

    cc_y = ey_y + 5 * ROW_PITCH + 34
    s.append(f'  <text x="{KX0 - 52}" y="{cc_y - 12}" font-size="{FS_GROUP}" '
             f'font-weight="bold" fill="#333333">Col-CEN (control)</text>')
    s += karyotype_group(cc_blocks, KBP_PER_PX["cc"], kb_per_px, KX0, cc_y,
                         cc_names)

    # scale bar, stated because the two source renders were drawn at
    # different kbp/px and had to be brought onto one scale
    sb_y = cc_y + 5 * ROW_PITCH + 20
    sb_w = 10000.0 / kb_per_px
    s.append(f'  <rect x="{KX0}" y="{sb_y}" width="{sb_w:.2f}" height="3" '
             f'fill="#444444"/>')
    s.append(f'  <text x="{KX0 + sb_w + 8:.2f}" y="{sb_y + 10}" font-size="{FS_BODY}" '
             f'fill="#555555">10 Mb, both genomes drawn at one scale</text>')

    # shared legend for the four labels the renderer drew
    lg_y = sb_y + 26
    lx = KX0
    for label in LEAF_ORDER:
        s.append(f'  <rect x="{lx}" y="{lg_y}" width="13" height="13" '
                 f'fill="{LEAF_COLOR[label]}" stroke="#666666" '
                 f'stroke-width="0.5"/>')
        s.append(f'  <text x="{lx + 19}" y="{lg_y + 11}" font-size="{FS_BODY}" '
                 f'fill="#555555">{esc(label)}</text>')
        lx += 19 + CHAR_W * len(label) + 24

    # ---- right: the composition the karyotype cannot show ------------------
    ey_r = comp[("ey15_2", "region")]
    cc_r = comp[("colcen", "region")]

    ry = top + 44
    s.append(f'  <text x="{CX0}" y="{ry}" font-size="{FS_BODY}" letter-spacing="0.8" '
             f'fill="#777777">REGION-SET COMPOSITION, ALL LABELS</text>')

    b1 = ry + 38
    s.append(f'  <text x="{CX0}" y="{b1 - 6}" font-size="{FS_BODY}" fill="#555555">'
             f'Ey15-2</text>')
    s += stacked_bar(ey_r, CX0, b1, CW, 26)
    b2 = b1 + 80
    s.append(f'  <text x="{CX0}" y="{b2 - 6}" font-size="{FS_BODY}" fill="#555555">'
             f'Col-CEN (control)</text>')
    s += stacked_bar(cc_r, CX0, b2, CW, 26)

    # bracket the non-leaf tail, which is too thin to read unaided
    ey_nl = nonleaf_total(ey_r)
    cc_nl = nonleaf_total(cc_r)
    nl_x = CX0 + CW * (100.0 - ey_nl) / 100.0
    s.append(f'  <path d="M {nl_x:.2f} {b1 + 30} L {nl_x:.2f} {b1 + 36} '
             f'L {CX0 + CW} {b1 + 36} L {CX0 + CW} {b1 + 30}" fill="none" '
             f'stroke="#B4741F" stroke-width="1.1"/>')
    s.append(f'  <text x="{CX0 + CW}" y="{b1 + 52}" font-size="{FS_BODY}" '
             f'fill="#B4741F" text-anchor="end">{ey_nl:.2f}% not resolved to a '
             f'leaf</text>')
    s.append(f'  <text x="{CX0 + CW}" y="{b2 + 46}" font-size="{FS_BODY}" '
             f'fill="#777777" text-anchor="end">{cc_nl:.3f}% in the control'
             f'</text>')

    ly = b2 + 66
    lx = CX0
    for label in LEAF_ORDER + NONLEAF_ORDER:
        fill = LEAF_COLOR.get(label) or NONLEAF_COLOR[label]
        item_w = 19 + CHAR_W * len(label) + 24
        if lx + item_w > CX0 + CW + 10:
            lx = CX0
            ly += 22
        s.append(f'  <rect x="{lx}" y="{ly}" width="13" height="13" '
                 f'fill="{fill}" stroke="#666666" stroke-width="0.5"/>')
        s.append(f'  <text x="{lx + 19}" y="{ly + 11}" font-size="{FS_BODY}" '
                 f'fill="#555555">{esc(label)}</text>')
        lx += item_w

    # ---- right, lower: where the unresolved mass actually goes -------------
    # Both destinations are drawn. The parent shift dominates in the chromosome
    # and region sets, but novel is the larger share in repeat and gene, so a
    # chart showing only the parent label would overstate the pattern.
    gy = ly + 56
    s.append(f'  <text x="{CX0}" y="{gy}" font-size="{FS_BODY}" letter-spacing="0.8" '
             f'fill="#777777">POSITIONS NOT RESOLVED TO A LEAF</text>')

    plot_y = gy + 20
    plot_h = 200
    axis_w = CW - 160
    ax0 = CX0 + 148
    vmax = 11.0
    for tick in (0, 5, 10):
        tx = ax0 + axis_w * tick / vmax
        s.append(f'  <line x1="{tx:.2f}" y1="{plot_y}" x2="{tx:.2f}" '
                 f'y2="{plot_y + plot_h}" stroke="#E2E2E2" stroke-width="1"/>')
        s.append(f'  <text x="{tx:.2f}" y="{plot_y + plot_h + 16}" '
                 f'font-size="{FS_SMALL}" fill="#888888" text-anchor="middle">'
                 f'{tick}%</text>')

    pitch = plot_h / len(SET_ORDER)
    for i, fset in enumerate(SET_ORDER):
        label = PARENT_LABEL[fset]
        ey_parent = comp[("ey15_2", fset)].get(label, 0.0)
        ey_novel = comp[("ey15_2", fset)].get("novel", 0.0)
        cc_parent = comp[("colcen", fset)].get(label, 0.0)
        row_y = plot_y + i * pitch + 2
        s.append(f'  <text x="{CX0}" y="{row_y + 13}" font-size="{FS_BODY}" '
                 f'fill="#555555">{esc(fset)}</text>')
        s.append(f'  <text x="{CX0}" y="{row_y + 28}" font-size="{FS_SMALL}" '
                 f'fill="#999999">{esc(label)}</text>')
        for j, (val, fill, tcol) in enumerate((
                (ey_parent, "#8A8A8A", "#666666"),
                (ey_novel, "#E8A33D", "#B4741F"),
                (cc_parent, "#C6C6C6", "#A0A0A0"))):
            by = row_y + j * 13
            bw = max(axis_w * val / vmax, 0.7)
            s.append(f'  <rect x="{ax0}" y="{by}" width="{bw:.2f}" height="10" '
                     f'fill="{fill}"/>')
            s.append(f'  <text x="{ax0 + bw + 6:.2f}" y="{by + 9.5}" '
                     f'font-size="{FS_SMALL}" fill="{tcol}">{fmt_pct(val)}</text>')

    key_y = plot_y + plot_h + 34
    kx = CX0
    for fill, text in (("#8A8A8A", "Ey15-2 parent"),
                       ("#E8A33D", "Ey15-2 novel"),
                       ("#C6C6C6", "Col-CEN parent")):
        s.append(f'  <rect x="{kx}" y="{key_y}" width="13" height="13" '
                 f'fill="{fill}"/>')
        s.append(f'  <text x="{kx + 19}" y="{key_y + 11}" font-size="{FS_BODY}" '
                 f'fill="#555555">{text}</text>')
        kx += 19 + CHAR_W * len(text) + 24

    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a3-dir", type=Path, default=None,
                    help="A3 cross-genome run directory (default: $KS_A3_DIR, "
                         "else <repo>/a3_crossgenome)")
    args = ap.parse_args()
    if args.a3_dir is not None:
        set_a3(args.a3_dir)
    for f in (EY_SVG, CC_SVG, COMP):
        if not f.is_file():
            raise SystemExit(f"panel C input missing: {f}")

    draft = DRAFT.read_text()
    body = draft.split("\n")
    # drop the closing tag and the opening svg/background, keep panels A and B
    start = next(i for i, ln in enumerate(body) if "Panel A" in ln)
    end = next(i for i, ln in enumerate(body) if ln.strip() == "</svg>")
    panels_ab = body[start:end]
    defs = "\n".join(body[1:start])

    comp = load_composition()
    panel_c = build_panel_c(comp)

    out = []
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">'
    )
    out.append(defs.replace(
        f'<rect x="0" y="0" width="1400" height="710" fill="#ffffff"/>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>'
    ))
    out += panels_ab
    out.append("")
    out.append("  <!-- ===================== Panel C: cross-accession transfer "
               "===================== -->")
    out += panel_c
    out.append("</svg>")

    OUT.write_text("\n".join(out) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
