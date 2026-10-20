#!/usr/bin/env python3
"""Write figures/fig7/build_workflow_panelsAB.svg, the hand-authored panels A and B
that scripts/make_build_figure.py reads and extends with the data-driven panel C.
Edit the constants and strings here, run this, then run make_build_figure.py.

Text sizes are SVG user units on the 1400-unit-wide canvas. The figure is
printed 183 mm wide, so 1 unit = 0.3705 pt; the FS_* constants keep body text
at >= 5 pt and panel letters at 8 pt at that size."""
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig7" / "build_workflow_panelsAB.svg"

Y0 = 400            # panel B offset
ROW = 151           # top of the build-step and database row
QY = 265            # top of the query assembly node
MONO = 'Menlo, Consolas, Courier New, monospace'; BLUE = "#1F5FBF"; SUB = "#444444"; HEAD = "#777777"
HX = ' letter-spacing="0.8"'
# Font sizes (user units; x 0.3705 = pt at 183 mm)
FS_LETTER = 21.6    # 8.0 pt, panel letters
FS_TITLE = 18.5     # 6.9 pt, panel titles
FS_NODE = 16        # 5.9 pt, node titles and hierarchy boxes
FS_EMPH = 15        # 5.6 pt, command names, resolution rules
FS_BODY = 14        # 5.2 pt, body text and headers
FS_SMALL = 13.5     # 5.0 pt, the annotate command label

def rect(x, y, w, h, fill="#ffffff", stroke="#555555", sw=1.4, rx=6):
    return f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
def t(x, y, s, size=FS_BODY, weight=None, fill="#222222", anchor="middle", family=None, extra=""):
    a = f' font-weight="{weight}"' if weight else ""; f = f' font-family="{family}"' if family else ""
    return f'  <text x="{x}" y="{y}" font-size="{size}"{a} fill="{fill}" text-anchor="{anchor}"{f}{extra}>{s}</text>'
def arrow(x1, y1, x2, y2):
    return f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#555555" stroke-width="1.6" marker-end="url(#arrow)"/>'
def node(x, y, w, h, title, lines, fill="#ffffff", stroke="#555555", sw=1.4, title_dy=None):
    n = len(lines)
    block = 22 + 16 * n                      # title line plus body lines
    top = y + (h - block) / 2
    out = [rect(x, y, w, h, fill, stroke, sw), t(x + w / 2, top + 16, title, FS_NODE, "bold")]
    for i, s in enumerate(lines):
        out.append(t(x + w / 2, top + 38 + 16 * i, s, fill=SUB))
    return out

L = ['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="710" viewBox="0 0 1400 710" font-family="Helvetica, Arial, sans-serif">',
     '  <defs>\n    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">\n      <path d="M 0 0 L 10 5 L 0 10 z" fill="#555555"/>\n    </marker>\n  </defs>\n',
     '  <rect x="0" y="0" width="1400" height="710" fill="#ffffff"/>\n',
     '  <!-- ===================== Panel A: build workflow ===================== -->']
L += [t(20, 36, "a", FS_LETTER, "bold", anchor="start"), t(50, 36, "Building a custom KaryoScope database", FS_TITLE, anchor="start")]
# left column
L += node(30, 86, 270, 100, "Raw annotation sources", ["RepeatMasker / EDTA tables,", "GFF3/GTF, cytobands,", "satellite monomer catalogues, .fai"])
L += [arrow(165, 186, 165, 236), t(177, 216, "prep-bed", FS_EMPH, fill=BLUE, anchor="start", family=MONO)]
L += node(30, 238, 270, 78, "Feature-set BED(s)", ["4th column = feature label,", "one BED per set"])
L += node(30, 330, 270, 52, "Genome FASTA", [])
L.append('  <path d="M312,238 L320,238 L320,382 L312,382" fill="none" stroke="#555555" stroke-width="1.4"/>')
cy = ROW + 50   # centre line of the build row
L.append(f'  <path d="M320,310 L338,310 L338,{cy} L358,{cy}" fill="none" stroke="#555555" stroke-width="1.6" marker-end="url(#arrow)"/>')
# build row
L += [f'  <path d="M362,{ROW-8} L362,{ROW-16} L982,{ROW-16} L982,{ROW-8}" fill="none" stroke="{BLUE}" stroke-width="1.4"/>',
      t(672, ROW - 24, "build", FS_EMPH, fill=BLUE, family=MONO)]
L += node(362, ROW, 170, 100, "Gap-fill", ["unannotated bases", "receive a background", "label"], fill="#F4F4F4")
L.append(arrow(532, cy, 570, cy))
L += node(572, ROW, 200, 100, "k-mer labels", ["a k-mer is recorded under", "every label where it occurs"], fill="#F4F4F4")
L.append(arrow(772, cy, 810, cy))
L += node(812, ROW, 170, 100, "Resolution", ["a k-mer with several", "labels resolves to the", "lowest common ancestor,", "or by priority"], fill="#F4F4F4")
L.append(arrow(982, cy, 1020, cy))
# right column
L += node(1022, ROW, 170, 100, "HKS database", ["base index", "one index per feature set", "manifest, hierarchy,", "colours"], stroke="#222222", sw=2)
L += node(1022, QY, 170, 52, "Query assembly", [])
bt, bb = ROW, QY + 52; bm = (bt + bb) / 2
L.append(f'  <path d="M1204,{bt} L1212,{bt} L1212,{bb} L1204,{bb}" fill="none" stroke="#555555" stroke-width="1.4"/>')
L += [arrow(1212, bm, 1279, bm), t(1244, bm - 9, "annotate", FS_SMALL, fill=BLUE, family=MONO)]
L += node(1282, bm - 36, 110, 72, "Annotations", ["base-pair", "resolution BEDs"])
# panel B
L += ['', '  <!-- ===================== Panel B: resolution example ===================== -->']
y = lambda v: v + Y0
L += [t(20, y(36), "b", FS_LETTER, "bold", anchor="start"), t(50, y(36), "Label resolution", FS_TITLE, anchor="start")]
L += [t(30, y(68), "HIERARCHY", fill=HEAD, anchor="start", extra=HX), t(420, y(68), "A K-MER WITH SEVERAL LABELS", fill=HEAD, anchor="start", extra=HX), t(990, y(68), "PRIORITY FILE (OPTIONAL)", fill=HEAD, anchor="start", extra=HX)]
L += [f'  <line x1="180" y1="{y(120)}" x2="180" y2="{y(162)}" stroke="#555555" stroke-width="1.4"/>',
      f'  <line x1="180" y1="{y(198)}" x2="115" y2="{y(232)}" stroke="#555555" stroke-width="1.4"/>',
      f'  <line x1="180" y1="{y(198)}" x2="245" y2="{y(232)}" stroke="#555555" stroke-width="1.4"/>']
for (x, yy, w, s) in [(120, 84, 120, "categorized"), (120, 162, 120, "Transposon"), (70, 232, 90, "LINE"), (200, 232, 90, "SINE")]:
    L += [rect(x, y(yy), w, 36), t(x + w / 2, y(yy + 23), s, FS_NODE, "bold")]
L.append(f'  <line x1="420" y1="{y(152)}" x2="960" y2="{y(152)}" stroke="#888888" stroke-width="1.2"/>')
L += [rect(440, y(136), 200, 32, fill="#F4F4F4", rx=4), rect(740, y(136), 200, 32, fill="#F4F4F4", rx=4)]
L += [f'  <rect x="580" y="{y(142)}" width="44" height="20" fill="#222222"/>', f'  <rect x="760" y="{y(142)}" width="44" height="20" fill="#222222"/>',
      f'  <path d="M602,{y(138)} C602,{y(98)} 782,{y(98)} 782,{y(138)}" fill="none" stroke="#222222" stroke-width="1.2" stroke-dasharray="3 3"/>',
      f'  <rect x="458" y="{y(142)}" width="44" height="20" fill="#ffffff" stroke="#222222" stroke-width="1.4"/>']
L += [t(692, y(92), "same k-mer"), t(540, y(192), "LINE, chr1:100&#8211;200", fill=SUB), t(840, y(192), "SINE, chr1:250&#8211;350", fill=SUB)]
L += [f'  <rect x="420" y="{y(222)}" width="24" height="13" fill="#ffffff" stroke="#222222" stroke-width="1.4"/>',
      t(456, y(233), 'found only in LINE &#8594; <tspan font-weight="bold">LINE</tspan>', FS_EMPH, anchor="start"),
      f'  <rect x="420" y="{y(246)}" width="24" height="13" fill="#222222"/>',
      t(456, y(257), 'found in LINE and SINE &#8594; <tspan font-weight="bold">Transposon</tspan>, their lowest common ancestor', FS_EMPH, anchor="start")]
L.append(rect(990, y(84), 380, 104, fill="#F4F4F4"))
for i, s in enumerate(["LINE        1  Transposon", "SINE        2  Transposon", "Transposon  1  categorized"]):
    L.append(t(1008, y(116 + 22 * i), s, FS_EMPH, anchor="start", family=MONO, extra=' xml:space="preserve"'))
L += [f'  <rect x="990" y="{y(206)}" width="24" height="13" fill="#222222"/>',
      t(1026, y(217), 'found in LINE and SINE &#8594; <tspan font-weight="bold">LINE</tspan>,', FS_EMPH, anchor="start"),
      t(1026, y(237), 'the higher-priority branch', FS_EMPH, anchor="start")]
L.append('</svg>')
OUT.write_text('\n'.join(L) + '\n')
print(f"wrote {OUT.relative_to(ROOT)}")
